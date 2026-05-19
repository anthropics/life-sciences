#!/usr/bin/env python3
"""build_igvreports.py — generic driver for the igv-reports skill.

Author: Samuel Ahuno
Purpose:
  Build self-contained HTML genomic-region reports with create_report
  (igv-reports). Two run modes:

    1. Single — direct CLI: --sites BED + --bam BAM(s) [+--vcf VCF]
       → one HTML at --output.

    2. Cohort — TSV samplesheet: one HTML per row + an index.html.
       Samplesheet columns (tab-separated, with header):
         sample   bam_tumor   bam_normal   vcf   sites_bed
       Optional fifth column: extra_tracks (comma-separated paths).

  Either way, the driver:
    - Resolves CpG islands, gencode, and RepeatMasker paths from
      databases_config.yaml for the chosen genome (skipping any not
      configured for that genome, with a warning).
    - Validates that the sites BED is headerless and well-formed.
    - Calls create_report with --flanking 300 --standalone by default.
    - Writes a logs/ entry capturing the resolved track list, the full
      command, the flanking value, and per-region embedded data sizes.

Usage:
  python build_igvreports.py --sites SITES.hg38.bed \\
      --bam tumor.bam normal.bam --vcf calls.vcf \\
      --genome hg38 --fasta /path/to/hg38.fa \\
      --no-default-tracks --output report.hg38.html

  python build_igvreports.py --samplesheet sheet.tsv \\
      --genome hg38 --fasta /path/to/hg38.fa \\
      --no-default-tracks --output-dir results/cohort/

Defaults can be shipped via a YAML at IGV_REPORTS_DB_CONFIG (see the
references/databases_config_paths.md schema), in which case --fasta /
--no-default-tracks are not needed.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

try:
    import yaml  # PyYAML
except ImportError:
    print("ERROR: PyYAML not available. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(2)

# Optional: point IGV_REPORTS_DB_CONFIG at a YAML mapping genome IDs to
# FASTA + CGI + gencode + RepeatMasker paths so --genome resolves tracks
# automatically. Without it, pass --fasta and --no-default-tracks explicitly.
_DB_CONFIG_ENV = os.environ.get("IGV_REPORTS_DB_CONFIG")
DEFAULT_DBCONFIG = Path(_DB_CONFIG_ENV) if _DB_CONFIG_ENV else None
DEFAULT_FLANKING = 300
# Optional: point IGV_REPORTS_SIF at an apptainer SIF for offline / HPC runs.
# Galaxy depot: https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7cba7a3_0
_SIF_ENV = os.environ.get("IGV_REPORTS_SIF")
IGVREPORTS_SIF = Path(_SIF_ENV) if _SIF_ENV else None


def apptainer_bind_args() -> list[str]:
    """Build `--bind <path>` tokens for singularity, skipping paths that don't
    exist. Source: `$IGV_REPORTS_BIND` (colon-separated). Empty / unset = no
    binds (singularity tolerates this; you only need binds when your data
    lives outside the container's default-visible filesystem)."""
    raw = os.environ.get("IGV_REPORTS_BIND", "")
    if not raw:
        return []
    tokens: list[str] = []
    for p in raw.split(":"):
        if p and Path(p).exists():
            tokens.extend(["--bind", p])
    return tokens

GENOME_ALIASES = {
    "hg38": "hg38",
    "GRCh38": "hg38",
    "mm10": "mm10",
    "GRCm38": "mm10",
    "mm39": "mm39",
    "GRCm39": "mm39",
    "t2t": "t2t_CHM13v2_plusY",
    "chm13": "t2t_CHM13v2_plusY",
    "T2T": "t2t_CHM13v2_plusY",
    "T2T-CHM13": "t2t_CHM13v2_plusY",
    "t2t_CHM13v2_plusY": "t2t_CHM13v2_plusY",
    "GRCh37": "GRCh37",
    "hg19": "GRCh37",
}


def setup_logger(log_path: Path) -> logging.Logger:
    """Dual-handler logger: file + stderr, with timestamp prefix."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    log = logging.getLogger("igv_reports")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    log.addHandler(fh)
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    log.addHandler(sh)
    return log


def resolve_genome(genome: str) -> str:
    canon = GENOME_ALIASES.get(genome)
    if not canon:
        raise SystemExit(
            f"ERROR: unknown genome '{genome}'. Supported: {sorted(set(GENOME_ALIASES.values()))}"
        )
    return canon


def load_db_config(path: Path) -> dict:
    """Load the databases YAML. Returns {} (with a warning to stderr) if the
    file is missing — callers must handle empty cfg gracefully.

    Users without a YAML can pass --fasta and --no-default-tracks
    on the driver, OR set $IGV_REPORTS_DB_CONFIG to their own YAML."""
    if not path.exists():
        sys.stderr.write(
            f"[build_igvreports] WARNING: db-config not found at {path}\n"
            "  Set $IGV_REPORTS_DB_CONFIG to point at your YAML, or pass\n"
            "  --fasta PATH and --no-default-tracks to bypass it entirely.\n"
        )
        return {}
    with path.open() as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg


def resolve_default_tracks(cfg: dict, genome: str, log: logging.Logger) -> list[str]:
    """Return ordered list of default tracks present on disk for this genome.

    Order matters — first entry renders at the bottom of the IGV.js view by
    default? Actually igv-reports renders --tracks in the order passed,
    top-to-bottom. We put annotation tracks LAST so they sit below the
    BAM/VCF data the user is actually inspecting.

    Empty cfg (e.g. generic, no databases YAML) → returns [] with a warning.
    """
    g = cfg.get("reference_genomes", {}).get("local", {}).get(genome, {})
    if not g:
        log.warning(
            f"no entry for genome '{genome}' in db-config — skipping default tracks. "
            "Pass --extra-track or --track-config for annotation tracks."
        )
        return []

    tracks: list[str] = []

    # CpG islands.
    cgi = g.get("CpGIslands")
    if cgi and Path(cgi).exists():
        tracks.append(cgi)
    else:
        log.warning(f"CpG islands track missing for {genome} (key=CpGIslands, value={cgi})")

    # Gencode. For hg38 prefer the bgzip+tabix .gff3.gz sibling if present.
    gtf = g.get("gtf")
    gencode_track: str | None = None
    if genome == "hg38" and gtf:
        sibling = (Path(gtf).parent / "gencode.v47.annotation.gff3.gz")
        if sibling.exists() and (sibling.parent / (sibling.name + ".tbi")).exists():
            gencode_track = str(sibling)
            log.info(f"  hg38: using full gencode annotation: {sibling}")
    if gencode_track is None and gtf and Path(gtf).exists():
        gencode_track = gtf
    if gencode_track:
        tracks.append(gencode_track)
    else:
        log.warning(f"Gencode track missing for {genome}")

    # RepeatMasker.
    rmsk = g.get("repMaskerBed")
    if rmsk and Path(rmsk).exists():
        tracks.append(rmsk)
    else:
        log.warning(f"RepeatMasker track not configured for {genome}")

    return tracks


def fasta_for(cfg: dict, genome: str) -> str:
    """Resolve a FASTA path from the db-config. Users without a
    YAML can bypass this by passing --fasta PATH on the driver."""
    try:
        fasta = cfg["reference_genomes"]["local"][genome].get("fasta")
    except (KeyError, TypeError):
        raise SystemExit(
            f"ERROR: db-config has no '{genome}' entry to resolve FASTA from.\n"
            "       Pass --fasta PATH explicitly, or set $IGV_REPORTS_DB_CONFIG\n"
            "       to a YAML that defines reference_genomes.local.<genome>.fasta."
        )
    if not fasta or not Path(fasta).exists():
        raise SystemExit(f"ERROR: FASTA missing for {genome}: {fasta}")
    if not Path(fasta + ".fai").exists():
        raise SystemExit(
            f"ERROR: FASTA index missing for {fasta} — run `samtools faidx {fasta}`"
        )
    return fasta


def validate_bams(bams: list[Path]) -> None:
    """Sanity-check BAMs have a coindex sibling (.bai or .csi).

    `create_report` needs a random-access index to slice BAMs at each region;
    a BAM with no sibling index produces an obscure pysam error several
    layers in. Catch it up front with an actionable message."""
    for bam in bams:
        if not bam.exists():
            raise SystemExit(f"ERROR: BAM not found: {bam}")
        if not (bam.with_suffix(bam.suffix + ".bai").exists()
                or bam.with_suffix(bam.suffix + ".csi").exists()
                or bam.with_suffix(".bai").exists()
                or bam.with_suffix(".csi").exists()):
            raise SystemExit(
                f"ERROR: BAM index missing for {bam} — create_report cannot slice it.\n"
                f"       Fix: samtools index {bam}\n"
                f"       (or `samtools index -c {bam}` for a .csi index on contigs >512 Mb)"
            )


def validate_sites_bed(bed: Path) -> None:
    """Sanity-check the sites BED before invoking create_report.

    create_report's BED parser is positional. It skips lines starting with
    `#` or `track ` (so the lab's `#chrom\\tstart\\tend\\tname` header is
    fine), but a non-comment header row like `chrom\\tstart\\tend` crashes
    with `ValueError: invalid literal for int()`. We mirror create_report's
    line-skipping logic and emit an informative error if any data row has
    non-numeric start/end."""
    if not bed.exists():
        raise SystemExit(f"ERROR: sites BED not found: {bed}")
    with bed.open() as fh:
        for i, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track "):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                raise SystemExit(f"ERROR: {bed}:{i}: BED needs >=3 tab-separated columns; got {cols!r}")
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                raise SystemExit(
                    f"ERROR: {bed}:{i}: non-numeric start/end — likely a header row.\n"
                    "       igv-reports' BED parser is positional and chokes on non-comment\n"
                    "       headers. Prefix the header with `#` (skipped by create_report\n"
                    "       and matches the lab's BED-output convention) or strip it."
                )
            if start >= end:
                raise SystemExit(f"ERROR: {bed}:{i}: start ({start}) >= end ({end})")


def find_create_report() -> str:
    """Resolve `create_report` on PATH (provided by `pip install igv-reports`
    or any conda env that activated it)."""
    cr = shutil.which("create_report")
    if cr:
        return cr
    raise SystemExit(
        "ERROR: create_report not on PATH.\n"
        "  Install:              pip install -U 'igv-reports>=1.16.0'\n"
        "  Offline / air-gapped: point IGV_REPORTS_SIF at an igv-reports SIF\n"
        "                        and rerun with --apptainer (Galaxy depot:\n"
        "                        https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7cba7a3_0)"
    )


def apptainer_create_report_prefix(sif: Path) -> list[str]:
    """Return the `singularity exec --cleanenv [--bind <path> ...] <sif>
    create_report` prefix. Used when --apptainer is passed; avoids the NFS
    conda cold-start tax (rules/apptainer_vs_conda.md). The default SIF is
    a dedicated igv-reports container (igv-reports_1.16.0.sif, ~83 MB)
    pulled from the Galaxy depot. Override via $IGV_REPORTS_SIF.

    --cleanenv: scrubs host env vars so they don't leak into the SIF.
    Specifically: host SSL_CERT_FILE / SSL_CERT_DIR on RHEL 8 point at paths
    that don't exist inside Galaxy-depot SIFs, and create_report's standalone-
    HTML build path performs an HTTPS GET (for the IGV.js ideogram or similar)
    that aborts with `[SSL: CERTIFICATE_VERIFY_FAILED]`. See
    rules/apptainer_env_leak.md.

    Binds: see `apptainer_bind_args()` — conditional on path existence."""
    if sif is None:
        raise SystemExit(
            "ERROR: --apptainer requested but $IGV_REPORTS_SIF is not set.\n"
            "       Set IGV_REPORTS_SIF to a SIF path and rerun, e.g.:\n"
            "         export IGV_REPORTS_SIF=/path/to/igv-reports_1.16.0.sif\n"
            "       Pull the SIF first if needed:\n"
            "         wget -O \"$IGV_REPORTS_SIF\" \\\n"
            "           'https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7cba7a3_0'"
        )
    if not sif.exists():
        raise SystemExit(
            f"ERROR: apptainer SIF not found: {sif}\n"
            "       Pull with one of:\n"
            f"         apptainer pull {sif} \\\n"
            "           docker://igv-org/igv-reports:1.16.0\n"
            f"         wget -O {sif} \\\n"
            "           'https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7e72e81_0'\n"
            "       Or set $IGV_REPORTS_SIF to a SIF you already have."
        )
    return ["singularity", "exec", "--cleanenv", *apptainer_bind_args(), str(sif), "create_report"]


def _read_sites_bed_rows(sites: Path) -> list[dict]:
    """Read a sites BED into a list of dicts, one per data row. Lines
    starting with `#` / `track ` / `browser ` are skipped (same rule
    `validate_sites_bed` uses). The 4th column (`name`) becomes the UID
    when present; otherwise an auto-generated `region_<idx>` is used."""
    rows: list[dict] = []
    with sites.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom, start_s, end_s = cols[0], cols[1], cols[2]
            name = cols[3].strip() if len(cols) >= 4 and cols[3].strip() else ""
            rows.append({
                "chrom": chrom,
                "start": int(start_s),
                "end": int(end_s),
                "name": name,
            })
    for idx, r in enumerate(rows, start=1):
        if not r["name"]:
            r["name"] = f"region_{idx:03d}"
        r["bed_row_idx"] = idx
    return rows


def _write_igver_regions_bed(rows: list[dict], flanking: int, out: Path) -> None:
    """Emit a BED with `--flanking` baked into start/end and UID in col 4.
    Filename collisions in igver's `chr-start-end.<uid>.png` are avoided by
    the auto-assigned UIDs in `_read_sites_bed_rows`."""
    with out.open("w") as fh:
        for r in rows:
            start = max(0, r["start"] - flanking)
            end = r["end"] + flanking
            fh.write(f"{r['chrom']}\t{start}\t{end}\t{r['name']}\n")


def _write_igver_input_list(tracks: list[str], out: Path) -> None:
    """One path per line — igver's `-i FOO.txt` consumes this verbatim."""
    with out.open("w") as fh:
        for t in tracks:
            fh.write(f"{t}\n")


def _resolve_igver_cmd(override: str | None) -> list[str]:
    """Return the argv prefix used to invoke igver. Resolution order:
    1. Explicit override (split on whitespace — supports `apptainer exec ... igver`).
    2. $IGVER_CMD env var (same shape as the override).
    3. `igver` on PATH."""
    if override:
        return override.split()
    env_cmd = os.environ.get("IGVER_CMD")
    if env_cmd:
        return env_cmd.split()
    on_path = shutil.which("igver")
    if on_path:
        return [on_path]
    raise SystemExit(
        "ERROR: igver not found.\n"
        "  Install: pip install igver\n"
        "  Override: --igver-cmd 'apptainer exec /path/to/igver.sif igver'\n"
        "  Or set $IGVER_CMD"
    )


def build_pngs_with_igver(
    sites: Path,
    tracks: list[str],
    genome: str,
    flanking: int,
    out_dir: Path,
    log: logging.Logger,
    html_path: Path,
    igver_cmd: str | None = None,
    dpi: int = 300,
    display_mode: str = "collapse",
    panel_height: int | None = None,
    fmt: str = "png",
) -> Path:
    """Invoke igver against the same sites + track list that drove
    create_report, write a manifest mapping each BED row to its PNG path
    and HTML row, return the manifest path.

    Consistency contract (the five levers from the design):
      1. Same sites BED + same `flanking` baked into the BED rows we pass.
      2. Same resolved track list.
      3. `display_mode` chosen to match HTML defaults (collapse).
      4. UID-based filenames let a user pair PNG ↔ HTML by string match.
      5. The manifest TSV is the audit trail; verify_cohort.py reads it.

    Output layout (caller controls `out_dir`):
      out_dir/
        igver_regions.bed     - flanked BED with UIDs in col 4 (igver -r)
        igver_input.txt       - track paths, one per line (igver -i)
        png/                  - actual PNGs (igver -o); filenames are
                                chr-start-end.<uid>.<png|svg|pdf>
        manifest.tsv          - cross-artifact bridge to the HTML
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    png_dir = out_dir / "png"
    png_dir.mkdir(parents=True, exist_ok=True)

    rows = _read_sites_bed_rows(sites)
    if not rows:
        raise SystemExit(f"ERROR: no data rows found in sites BED: {sites}")

    regions_bed = out_dir / "igver_regions.bed"
    _write_igver_regions_bed(rows, flanking, regions_bed)
    input_txt = out_dir / "igver_input.txt"
    _write_igver_input_list(tracks, input_txt)

    cmd = list(_resolve_igver_cmd(igver_cmd)) + [
        "-i", str(input_txt),
        "-r", str(regions_bed),
        "-o", str(png_dir),
        "-g", genome,
        "-d", display_mode,
        "--dpi", str(dpi),
        "-f", fmt,
        "--no-singularity",
    ]
    if panel_height is not None:
        cmd.extend(["-p", str(panel_height)])

    log.info(f"  igver: dpi={dpi} display={display_mode} fmt={fmt} regions={len(rows)}")
    log.info(f"  igver cmd: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(f"igver FAILED for {sites}")
        log.error(f"stdout: {proc.stdout}")
        log.error(f"stderr: {proc.stderr}")
        raise SystemExit(proc.returncode)

    # PNG filename convention is set by igver's _parse_bed_file:
    # `<chrom>-<start>-<end>.<uid>.<ext>`. We reconstruct it here.
    ext = "svg" if fmt in ("svg", "pdf") else fmt
    manifest = out_dir / "manifest.tsv"
    with manifest.open("w") as fh:
        fh.write(
            "#bed_row_idx\tuid\tchrom\tstart_orig\tend_orig\t"
            "start_flanked\tend_flanked\tregion\tpng_path\thtml_path\thtml_table_row\n"
        )
        for r in rows:
            start_f = max(0, r["start"] - flanking)
            end_f = r["end"] + flanking
            fname = f"{r['chrom']}-{start_f}-{end_f}.{r['name']}.{ext}"
            png_rel = (png_dir / fname).resolve()
            html_rel = html_path.resolve()
            fh.write(
                f"{r['bed_row_idx']}\t{r['name']}\t{r['chrom']}\t"
                f"{r['start']}\t{r['end']}\t{start_f}\t{end_f}\t"
                f"{r['chrom']}:{start_f}-{end_f}\t{png_rel}\t{html_rel}\t"
                f"{r['bed_row_idx']}\n"
            )

    log.info(f"  png manifest: {manifest} ({len(rows)} rows)")

    # Inline existence check — igver exits 0 even when it fails to render
    # (silent exit-0 failure documented in the upstream skill's notes).
    # We can't trust the exit code, so verify every expected PNG path is
    # on disk and non-empty before returning. Without this a "successful"
    # build silently ships an empty png/ dir.
    missing: list[str] = []
    empty: list[str] = []
    for r in rows:
        start_f = max(0, r["start"] - flanking)
        end_f = r["end"] + flanking
        fname = f"{r['chrom']}-{start_f}-{end_f}.{r['name']}.{ext}"
        p = png_dir / fname
        if not p.exists():
            missing.append(fname)
        elif p.stat().st_size == 0:
            empty.append(fname)
    if missing or empty:
        log.error(
            f"igver returned exit 0 but {len(missing)} expected PNG(s) are missing "
            f"and {len(empty)} are zero-byte (out of {len(rows)} regions). "
            "This is a documented silent-failure mode of igver."
        )
        if missing:
            log.error(f"  missing: {missing[:5]}{'...' if len(missing) > 5 else ''}")
        if empty:
            log.error(f"  empty:   {empty[:5]}{'...' if len(empty) > 5 else ''}")
        raise SystemExit(
            f"ERROR: igver produced {len(rows) - len(missing) - len(empty)} of "
            f"{len(rows)} PNGs (silent exit-0 failure). Check the igver install "
            "path — `pip install igver` egg-link lacks the IGV Java binary; use "
            "an apptainer SIF via --igver-cmd or $IGVER_CMD."
        )

    return manifest


def build_one(
    sites: Path,
    bams: list[Path],
    vcf: Path | None,
    extra_tracks: list[Path],
    fasta: str,
    default_tracks: list[str],
    output: Path,
    title: str,
    flanking: int,
    log: logging.Logger,
    track_config: Path | None = None,
    report_type: str | None = None,
    info_columns: list[str] | None = None,
    use_apptainer: bool = False,
    also_png: bool = False,
    igver_cmd: str | None = None,
    png_dpi: int = 300,
    png_display_mode: str = "collapse",
    png_out_dir: Path | None = None,
) -> Path:
    """Run create_report for one site set and return the HTML path.

    Two track modes:
      * Default — positional `--tracks <path> <path> ...`. Used when
        `track_config` is None. BAM + VCF + extra + default annotations,
        in render order top-to-bottom.
      * track-config — `--track-config <json>`. Used when `track_config`
        is provided. The JSON is the source of truth; default_tracks,
        bams, vcf, extra_tracks are IGNORED (they go in the JSON instead).
        This is the path required for ONT methylation viewers (named
        tracks, per-track color/min/max/colorBy/displayMode).
    """
    validate_sites_bed(sites)
    # Only validate BAMs on the positional --tracks path. The --track-config
    # JSON has its own track-resolution semantics and may reference BAMs by
    # arbitrary url:; create_report itself will fail loudly there if needed.
    if track_config is None:
        validate_bams(bams)
    output.parent.mkdir(parents=True, exist_ok=True)

    create_report_cmd = (
        apptainer_create_report_prefix(IGVREPORTS_SIF) if use_apptainer
        else [find_create_report()]
    )

    cmd: list[str] = list(create_report_cmd) + [
        str(sites),
        "--fasta", fasta,
        "--flanking", str(flanking),
    ]

    if track_config is not None:
        cmd.extend(["--track-config", str(track_config)])
        log.info(f"  track-config: {track_config}  (defaults+bams+vcf bypassed)")
        if bams or vcf or extra_tracks or default_tracks:
            log.warning(
                "--track-config supplied; ignoring --bam/--vcf/--extra-track and "
                "auto-resolved default tracks. Put everything in the JSON instead."
            )
    else:
        # Track ordering: BAMs (data) -> VCF (calls) -> extra -> defaults (annotation, last).
        tracks: list[str] = [str(b) for b in bams]
        if vcf:
            tracks.append(str(vcf))
        tracks.extend(str(t) for t in extra_tracks)
        tracks.extend(default_tracks)
        cmd.extend(["--tracks", *tracks])
        log.info(f"  tracks (in render order):")
        for i, t in enumerate(tracks, start=1):
            log.info(f"    {i:>2}. {t}")

    if report_type:
        cmd.extend(["--type", report_type])
    if info_columns:
        cmd.extend(["--info-columns", *info_columns])

    cmd.extend([
        "--standalone",
        "--title", title,
        "--output", str(output),
    ])

    log.info(f"  cmd: {' '.join(cmd)}")
    log.info(f"  flanking_bp: {flanking}")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log.error(f"create_report FAILED for {sites}")
        log.error(f"stdout: {proc.stdout}")
        log.error(f"stderr: {proc.stderr}")
        raise SystemExit(proc.returncode)

    if output.exists():
        log.info(f"  HTML: {output} ({output.stat().st_size / 1024 / 1024:.2f} MB)")

    # PNG sidecar — same regions, same tracks, written next to the HTML.
    # On the --track-config path we extract every `url` from the JSON
    # (file resources only — http(s) URLs are skipped since igver can't
    # consume them); on the positional path we reuse the same ordered
    # list we just passed to create_report.
    if also_png:
        if track_config is not None:
            try:
                with track_config.open() as fh:
                    cfg_tracks = json.load(fh)
            except Exception as e:
                log.warning(f"--also-png: unable to parse --track-config JSON: {e} — skipping PNG step")
                return output
            png_tracks: list[str] = []
            for t in cfg_tracks if isinstance(cfg_tracks, list) else []:
                url = t.get("url") if isinstance(t, dict) else None
                if url and not str(url).startswith(("http://", "https://")):
                    png_tracks.append(str(url))
            if not png_tracks:
                log.warning("--also-png: track-config has no local-path tracks — skipping PNG step")
                return output
        else:
            png_tracks = [str(b) for b in bams]
            if vcf:
                png_tracks.append(str(vcf))
            png_tracks.extend(str(t) for t in extra_tracks)
            png_tracks.extend(default_tracks)

        out_dir = png_out_dir if png_out_dir is not None else (
            output.parent / f"png_{output.stem}"
        )
        parts = output.stem.split(".")
        genome_tag = parts[-1] if len(parts) >= 2 else "hg38"
        build_pngs_with_igver(
            sites=sites,
            tracks=png_tracks,
            genome=genome_tag,
            flanking=flanking,
            out_dir=out_dir,
            log=log,
            html_path=output,
            igver_cmd=igver_cmd,
            dpi=png_dpi,
            display_mode=png_display_mode,
        )

    return output


def parse_samplesheet(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open() as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for ln in fh:
            cols = ln.rstrip("\n").split("\t")
            if not cols or not cols[0].strip():
                continue
            row = dict(zip(header, cols))
            rows.append(row)
    required = {"sample", "sites_bed"}
    if rows and not required.issubset(rows[0].keys()):
        raise SystemExit(
            f"ERROR: samplesheet must have columns: sample, sites_bed (got {list(rows[0].keys())}).\n"
            "       Optional columns: bam_tumor, bam_normal, vcf, extra_tracks (comma-separated)."
        )
    return rows


def derive_log_path(out_dir: Path, override: Path | None = None) -> Path:
    """Choose a log dir matching the lab's `results/<run>/{reports,logs}/`
    sibling layout when possible. Fall back to `out_dir/logs/` (in-dir) when
    the sibling can't be created — `out_dir.parent` is root, read-only, or
    otherwise unwritable. Honor an explicit `override` unconditionally."""
    if override is not None:
        log_dir = override
    else:
        out_dir = out_dir.resolve()
        sibling = out_dir.parent / "logs"
        try:
            sibling.mkdir(parents=True, exist_ok=True)
            log_dir = sibling
        except (PermissionError, OSError):
            log_dir = out_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir / f"run_{datetime.now():%Y%m%d_%H%M%S}.log"


def write_index(report_paths: dict[str, Path], out: Path, title: str) -> Path:
    items = "\n".join(
        f'  <li><a href="{p.name}">{s}</a></li>'
        for s, p in sorted(report_paths.items())
    )
    out.write_text(
        "<!doctype html>\n<html><head><title>"
        + title
        + "</title><style>body{font-family:Arial,sans-serif;margin:2em}li{margin:0.4em 0}</style></head><body>\n"
        f"<h1>{title}</h1>\n<ul>\n{items}\n</ul>\n</body></html>\n"
    )
    return out


def run_anchors_generate(
    samplesheet: Path,
    sites_files: list[Path],
    out: Path,
    fail_on_fail: bool,
    log: logging.Logger,
) -> None:
    """Invoke `verify_anchors.py generate` once per distinct sites BED in the
    cohort, merging into a single anchors TSV at `out`. Most cohorts share
    one sites BED so this collapses to a single call; multi-sites cohorts
    get one anchor block per sites file."""
    script = Path(__file__).resolve().parent / "verify_anchors.py"
    if not script.exists():
        log.warning(f"anchors generate: script not found at {script} — skipping")
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("")  # truncate; per-sites blocks appended below
    for i, sites in enumerate(sites_files):
        block = out.with_suffix(f".part{i}.tsv")
        cmd = [
            sys.executable, str(script), "generate",
            "--samplesheet", str(samplesheet),
            "--sites", str(sites),
            "--out", str(block),
        ]
        log.info(f"anchors generate: {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        for line in (proc.stdout or "").splitlines():
            log.info(f"  anchors > {line}")
        if proc.stderr:
            for line in proc.stderr.splitlines():
                log.info(f"  anchors (stderr) > {line}")
        if proc.returncode != 0:
            if fail_on_fail:
                raise SystemExit(proc.returncode)
            log.warning(f"anchors generate exited {proc.returncode}; continuing")
            continue
        # Merge: keep header from first block, body rows from all.
        if i == 0:
            out.write_text(block.read_text())
        else:
            with out.open("a") as fh:
                for ln in block.read_text().splitlines():
                    if not ln or ln.startswith("#"):
                        continue
                    fh.write(ln + "\n")
        block.unlink()
    log.info(f"anchors generate: wrote {out}")


def run_anchors_verify(
    samplesheet: Path,
    reports_dir: Path,
    genome: str,
    anchors: Path,
    fail_on_fail: bool,
    log: logging.Logger,
) -> None:
    """Invoke `verify_anchors.py verify-cohort` after a cohort build."""
    script = Path(__file__).resolve().parent / "verify_anchors.py"
    if not script.exists():
        log.warning(f"anchors verify: script not found at {script} — skipping")
        return
    if not anchors.exists():
        log.warning(f"anchors verify: anchors TSV missing: {anchors} — skipping")
        return
    out = reports_dir / "cohort_verify_anchors.tsv"
    cmd = [
        sys.executable, str(script), "verify-cohort",
        "--samplesheet", str(samplesheet),
        "--reports-dir", str(reports_dir),
        "--genome", genome,
        "--anchors", str(anchors),
        "--out", str(out),
    ]
    if fail_on_fail:
        cmd.append("--fail-on-fail")
    log.info(f"anchors verify: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    for line in (proc.stdout or "").splitlines():
        log.info(f"  anchors > {line}")
    if proc.stderr:
        for line in proc.stderr.splitlines():
            log.info(f"  anchors (stderr) > {line}")
    log.info(f"anchors verify: TSV={out} exit={proc.returncode}")
    if proc.returncode != 0:
        if fail_on_fail:
            raise SystemExit(proc.returncode)
        log.warning(f"anchors verify exited {proc.returncode}; --fail-on-fail not set, continuing")


def run_cohort_verify(
    samplesheet: Path,
    reports_dir: Path,
    genome: str,
    db_config: Path,
    fail_on_fail: bool,
    log: logging.Logger,
) -> None:
    """Invoke verify_cohort.py at the end of a cohort build. Writes the TSV +
    summary next to the cohort's index.html. Fails the build if
    `fail_on_fail` is set and the verifier exits nonzero."""
    verify_script = Path(__file__).resolve().parent / "verify_cohort.py"
    if not verify_script.exists():
        log.warning(f"verify_cohort: script not found at {verify_script} — skipping")
        return
    tsv_out = reports_dir / "cohort_verify.tsv"
    md_out = reports_dir / "cohort_verify.summary.md"
    cmd = [
        sys.executable, str(verify_script),
        "--samplesheet", str(samplesheet),
        "--reports-dir", str(reports_dir),
        "--genome", genome,
        "--db-config", str(db_config),
        "--out", str(tsv_out),
        "--summary", str(md_out),
    ]
    if fail_on_fail:
        cmd.append("--fail-on-fail")
    log.info(f"verify_cohort: running {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    # Mirror the verifier's stdout/stderr into the run log so audit-trail stays single-source.
    for line in (proc.stdout or "").splitlines():
        log.info(f"  verify_cohort > {line}")
    if proc.stderr:
        for line in proc.stderr.splitlines():
            log.warning(f"  verify_cohort (stderr) > {line}")
    log.info(f"verify_cohort: TSV={tsv_out} summary={md_out} exit={proc.returncode}")
    if proc.returncode != 0:
        if fail_on_fail:
            raise SystemExit(proc.returncode)
        log.warning(f"verify_cohort: exited {proc.returncode} but --fail-on-fail not set; continuing")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="hg38 | mm10 | mm39 | t2t | GRCh37 (alias-tolerant)")
    ap.add_argument("--db-config", default=(str(DEFAULT_DBCONFIG) if DEFAULT_DBCONFIG else None), help=(
        "YAML resolving genome -> {fasta, CpGIslands, gtf, repMaskerBed}. "
        "Schema: see references/databases_config_paths.md. "
        "Override via $IGV_REPORTS_DB_CONFIG, or skip entirely with --fasta + --no-default-tracks."
    ))
    ap.add_argument("--fasta", help=(
        "Explicit FASTA path; bypasses --db-config for FASTA lookup. "
        "Required when --db-config is not set or doesn't list the chosen genome. "
        "Requires a sibling .fai (run `samtools faidx`)."
    ))
    ap.add_argument("--no-default-tracks", action="store_true", help=(
        "Skip the CpG-islands/gencode/RepeatMasker auto-tracks from --db-config. "
        "Combine with --fasta and --extra-track to operate without a databases YAML."
    ))
    ap.add_argument("--flanking", type=int, default=DEFAULT_FLANKING)
    ap.add_argument("--extra-track", action="append", default=[], help="(repeat) extra track path; rendered above default annotations")

    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--samplesheet", help="TSV: sample, [bam_tumor, bam_normal, vcf,] sites_bed[, extra_tracks]")
    mode.add_argument("--sites", help="path to sites BED for single-sample mode")

    ap.add_argument("--bam", nargs="*", default=[], help="BAM/CRAM tracks (single-sample mode)")
    ap.add_argument("--vcf", help="VCF track (single-sample mode)")

    ap.add_argument("--output", help="output HTML path (single-sample mode)")
    ap.add_argument("--output-dir", help="output dir for cohort mode (default: ./reports)")
    ap.add_argument("--title", default=None, help="report title; defaults to sample name + genome")

    ap.add_argument(
        "--track-config",
        help="path to a tracks.json (igv.js track config). When set, the JSON is "
             "passed straight to create_report --track-config and all default "
             "tracks / --bam / --vcf / --extra-track are bypassed. Use this for "
             "ONT methylation viewers — see examples/methylation_ont/.",
    )
    ap.add_argument(
        "--type",
        dest="report_type",
        choices=["mutation", "fusion", "junction"],
        default=None,
        help="create_report --type. Sets viewer behaviour at each site.",
    )
    ap.add_argument(
        "--info-columns",
        nargs="*",
        default=[],
        help="VCF INFO or BED columns to surface in the variant table. "
             "For BED sites, 'name' is the most useful.",
    )
    ap.add_argument(
        "--also-png",
        action="store_true",
        help="After create_report finishes, invoke igver against the same "
             "sites BED + track list to produce per-region PNGs alongside "
             "the HTML. PNGs land in <html_parent>/png_<sample>/png/ with "
             "filename `<chr-start-end>.<uid>.png` (uid = BED `name` col, "
             "auto-assigned `region_<idx>` when missing). A manifest TSV "
             "bridges PNG ↔ HTML rows. Requires `igver` on PATH or "
             "$IGVER_CMD / --igver-cmd override.",
    )
    ap.add_argument(
        "--igver-cmd",
        default=None,
        help="Override the igver invocation. Resolution order: this flag, "
             "$IGVER_CMD, `igver` on PATH. Pass the full command including "
             "any apptainer wrapper, e.g. 'apptainer exec /path/to/igver.sif igver'.",
    )
    ap.add_argument(
        "--png-dpi",
        type=int,
        default=300,
        help="DPI for igver PNG output (default 300; bump to 600 for slide-quality).",
    )
    ap.add_argument(
        "--png-display-mode",
        choices=["expand", "collapse", "squish"],
        default="collapse",
        help="igver `-d` flag. Default 'collapse' to match the HTML's BAM "
             "BAM_DEFAULTS displayMode. Use 'expand' for per-read SV inspection.",
    )
    ap.add_argument(
        "--apptainer",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Run create_report from inside the apptainer SIF pointed to "
             "by $IGV_REPORTS_SIF (dedicated igv-reports 1.16.0 SIF, ~83 MB; "
             "pull from the Galaxy depot). Skips the NFS conda cold-start "
             "tax on HPC. Default: auto-detect — on if SLURM_JOB_ID is set "
             "AND $IGV_REPORTS_SIF points at an existing SIF, off otherwise. "
             "Override either way with --apptainer / --no-apptainer.",
    )
    ap.add_argument(
        "--log-dir",
        help="explicit log directory. Default: sibling 'logs/' of the output "
             "dir (matches results/<run>/{reports,logs}/ lab layout); falls "
             "back to <out_dir>/logs/ when the sibling is unwritable.",
    )
    ap.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        help="Number of parallel per-sample builds in cohort (--samplesheet) "
             "mode. Each worker invokes create_report in a subprocess, so the "
             "win comes from running multiple slicers concurrently against "
             "different BAMs. I/O-bound on the BAM-slice step, so threads "
             "scale well to ~min(N_samples, N_cores). Default 1 (sequential, "
             "preserves prior behavior). Has no effect in single-sample mode.",
    )
    ap.add_argument(
        "--verify",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run scripts/verify_cohort.py at the end of cohort builds "
             "(--samplesheet mode). Single-sample (--sites) mode is unaffected "
             "and emits no cohort verify TSV. Default: on. Use --no-verify to "
             "skip. The verifier inherits --fail-on-fail.",
    )
    ap.add_argument(
        "--fail-on-fail",
        action="store_true",
        help="Propagated to verify_cohort.py and verify_anchors.py: exit "
             "nonzero if any verifier check is FAIL. Only meaningful with "
             "--verify / --anchors-mode and --samplesheet.",
    )
    ap.add_argument(
        "--anchors-mode",
        choices=["off", "generate", "verify"],
        default="off",
        help="Content (read-count) verification — opt-in because it shells "
             "out to samtools per (sample, region) and is slow. 'generate' "
             "runs samtools view -c against source BAMs at build time and "
             "freezes the counts to --anchors (becomes regression fixture). "
             "'verify' decodes each BAM slice from the built HTMLs and "
             "compares to --anchors. 'off' (default) skips. See "
             "examples/anchor_verify_demo/.",
    )
    ap.add_argument(
        "--anchors",
        help="Path to anchors TSV. With --anchors-mode generate: output. "
             "With --anchors-mode verify: input. Ignored when mode=off.",
    )

    args = ap.parse_args()

    genome = resolve_genome(args.genome)
    # Only load db-config when something actually needs it (fasta lookup or
    # default tracks). Saves the warning noise + lets a fully-explicit
    # --fasta + --no-default-tracks invocation run with no YAML at all.
    need_db_config = (not args.fasta) or (not args.no_default_tracks)
    if need_db_config and not args.db_config:
        raise SystemExit(
            "ERROR: need a databases YAML to resolve FASTA / default tracks.\n"
            "  Either:\n"
            "    - pass --fasta /path/to/genome.fa --no-default-tracks "
            "[--extra-track ...], or\n"
            "    - set $IGV_REPORTS_DB_CONFIG / --db-config to a YAML matching\n"
            "      the schema in references/databases_config_paths.md."
        )
    cfg = load_db_config(Path(args.db_config)) if need_db_config else {}
    if args.fasta:
        fasta = args.fasta
        if not Path(fasta).exists():
            raise SystemExit(f"ERROR: --fasta path not found: {fasta}")
        if not Path(fasta + ".fai").exists():
            raise SystemExit(
                f"ERROR: FASTA index missing for {fasta} — run `samtools faidx {fasta}`"
            )
    else:
        fasta = fasta_for(cfg, genome)

    # Logger placed alongside the output. See derive_log_path docstring.
    if args.samplesheet:
        out_dir = Path(args.output_dir or "reports")
    else:
        if not args.output:
            raise SystemExit("ERROR: --output required in single-sample mode")
        out_dir = Path(args.output).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = derive_log_path(out_dir, Path(args.log_dir) if args.log_dir else None)
    log = setup_logger(log_path)

    log.info(f"=== igv-reports skill, genome={genome} ===")
    log.info(f"db_config: {args.db_config}")
    log.info(f"fasta:     {fasta}")
    log.info(f"flanking:  {args.flanking} bp (default {DEFAULT_FLANKING})")

    # Resolve --apptainer auto-detect. Tri-state:
    #   user said --apptainer        -> True
    #   user said --no-apptainer     -> False
    #   user said nothing (None)     -> True iff SLURM_JOB_ID is in env
    # Rationale: on a fresh SLURM compute node, the NFS conda cold-start tax
    # (~1-2 M page faults, ~2.5 us each) is large; the dedicated SIF skips it.
    # On the login node, conda is usually warm and the simpler path wins.
    # See rules/apptainer_vs_conda.md.
    slurm_job = os.environ.get("SLURM_JOB_ID")
    if args.apptainer is None:
        # Auto-enable SIF mode only when both (a) we're on a SLURM compute
        # node where the conda cold-start tax bites, AND (b) $IGV_REPORTS_SIF
        # points at an existing SIF. The existence check protects users from
        # a confusing SIF-not-found error when they didn't ask for apptainer.
        sif_ok = IGVREPORTS_SIF is not None and IGVREPORTS_SIF.exists()
        args.apptainer = bool(slurm_job) and sif_ok
        if args.apptainer:
            decision = f"auto-enabled (SLURM_JOB_ID={slurm_job}, SIF={IGVREPORTS_SIF})"
        elif slurm_job and IGVREPORTS_SIF is None:
            decision = (
                f"auto-disabled (SLURM_JOB_ID={slurm_job} set, but $IGV_REPORTS_SIF unset; "
                f"falling back to PATH create_report)"
            )
        elif slurm_job:
            decision = (
                f"auto-disabled (SLURM_JOB_ID={slurm_job} set, but SIF not found at "
                f"{IGVREPORTS_SIF}; falling back to PATH create_report)"
            )
        else:
            decision = "auto-disabled (no SLURM_JOB_ID; PATH create_report path)"
        log.info(f"apptainer: {decision}")
    else:
        log.info(f"apptainer: {args.apptainer} (explicit)")

    if args.no_default_tracks:
        default_tracks: list[str] = []
        log.info("default tracks: skipped (--no-default-tracks)")
    else:
        default_tracks = resolve_default_tracks(cfg, genome, log)
        log.info(f"default tracks resolved: {len(default_tracks)}")
        for t in default_tracks:
            log.info(f"  - {t}")

    extra_tracks = [Path(p) for p in args.extra_track]

    track_config = Path(args.track_config) if args.track_config else None
    if track_config is not None and not track_config.exists():
        raise SystemExit(f"ERROR: --track-config file not found: {track_config}")

    if args.sites:
        title = args.title or f"{Path(args.sites).stem} ({genome})"
        build_one(
            sites=Path(args.sites),
            bams=[Path(b) for b in args.bam],
            vcf=Path(args.vcf) if args.vcf else None,
            extra_tracks=extra_tracks,
            fasta=fasta,
            default_tracks=default_tracks,
            output=Path(args.output),
            title=title,
            flanking=args.flanking,
            log=log,
            track_config=track_config,
            report_type=args.report_type,
            info_columns=args.info_columns,
            use_apptainer=args.apptainer,
            also_png=args.also_png,
            igver_cmd=args.igver_cmd,
            png_dpi=args.png_dpi,
            png_display_mode=args.png_display_mode,
        )
    else:
        rows = parse_samplesheet(Path(args.samplesheet))
        n_jobs = max(1, args.jobs)
        n_workers = min(n_jobs, len(rows)) if rows else 1
        mode = "sequential" if n_workers == 1 else f"parallel ({n_workers} workers)"
        log.info(f"cohort: {len(rows)} samples from {args.samplesheet} — {mode}")
        report_paths: dict[str, Path] = {}
        failures: list[tuple[str, str]] = []  # (sample, error_message)

        def _build_row(row: dict) -> tuple[str, Path]:
            """Build one sample. Runs in a worker thread when --jobs > 1.

            Returns (sample, out_html). Raises on build failure — caught by
            the executor and surfaced via future.exception() in the caller."""
            sample = row["sample"]
            sites = Path(row["sites_bed"])
            bams = [Path(row[k]) for k in ("bam_tumor", "bam_normal") if row.get(k)]
            vcf = Path(row["vcf"]) if row.get("vcf") else None
            sample_extras = list(extra_tracks)
            if row.get("extra_tracks"):
                sample_extras += [Path(p.strip()) for p in row["extra_tracks"].split(",") if p.strip()]
            out_html = out_dir / f"{sample}.{genome}.html"
            title = args.title or f"{sample} ({genome})"
            log.info(f"=== {sample} ===")
            build_one(
                sites=sites, bams=bams, vcf=vcf, extra_tracks=sample_extras,
                fasta=fasta, default_tracks=default_tracks,
                output=out_html, title=title, flanking=args.flanking, log=log,
                track_config=track_config,
                report_type=args.report_type,
                info_columns=args.info_columns,
                use_apptainer=args.apptainer,
                also_png=args.also_png,
                igver_cmd=args.igver_cmd,
                png_dpi=args.png_dpi,
                png_display_mode=args.png_display_mode,
            )
            return sample, out_html

        # ThreadPoolExecutor is the right primitive here: build_one() spends
        # nearly all its wall time inside subprocess.run(create_report), which
        # releases the GIL — so threads scale linearly to the number of
        # concurrent create_report processes the host can support. Don't use
        # ProcessPoolExecutor: build_one() captures a non-picklable logger.
        if n_workers == 1:
            for row in rows:
                try:
                    sample, out_html = _build_row(row)
                    report_paths[sample] = out_html
                except SystemExit as exc:
                    failures.append((row.get("sample", "?"), f"exit={exc.code}"))
                except Exception as exc:
                    failures.append((row.get("sample", "?"), f"{type(exc).__name__}: {exc}"))
        else:
            with ThreadPoolExecutor(max_workers=n_workers) as pool:
                future_to_sample = {
                    pool.submit(_build_row, row): row.get("sample", "?") for row in rows
                }
                # as_completed lets failures surface immediately while other
                # samples continue building. We collect all errors and decide
                # whether to fail the whole run at the end.
                for fut in as_completed(future_to_sample):
                    sample_name = future_to_sample[fut]
                    try:
                        sample, out_html = fut.result()
                        report_paths[sample] = out_html
                    except SystemExit as exc:
                        failures.append((sample_name, f"exit={exc.code}"))
                    except Exception as exc:
                        failures.append((sample_name, f"{type(exc).__name__}: {exc}"))

        if failures:
            log.error(f"cohort: {len(failures)} of {len(rows)} samples FAILED:")
            for s, err in failures:
                log.error(f"  - {s}: {err}")
            # Always raise on build failures — these aren't verifier soft-fails,
            # they're missing HTMLs. --fail-on-fail is for verifier behavior.
            raise SystemExit(1)

        idx = write_index(report_paths, out_dir / "index.html", f"igv-reports cohort ({genome})")
        log.info(f"Wrote cohort index: {idx}")

        if args.verify:
            run_cohort_verify(
                samplesheet=Path(args.samplesheet),
                reports_dir=out_dir,
                genome=genome,
                db_config=Path(args.db_config),
                fail_on_fail=args.fail_on_fail,
                log=log,
            )
        else:
            log.info("verify_cohort: skipped (--no-verify)")

        if args.anchors_mode != "off":
            if not args.anchors:
                raise SystemExit("ERROR: --anchors PATH required when --anchors-mode != off")
            anchors_path = Path(args.anchors)
            if args.anchors_mode == "generate":
                sites_files = sorted({Path(r["sites_bed"]) for r in rows if r.get("sites_bed")})
                run_anchors_generate(
                    samplesheet=Path(args.samplesheet),
                    sites_files=sites_files,
                    out=anchors_path,
                    fail_on_fail=args.fail_on_fail,
                    log=log,
                )
            else:  # verify
                run_anchors_verify(
                    samplesheet=Path(args.samplesheet),
                    reports_dir=out_dir,
                    genome=genome,
                    anchors=anchors_path,
                    fail_on_fail=args.fail_on_fail,
                    log=log,
                )

    log.info(f"=== DONE: build_igvreports.py completed successfully ===")


if __name__ == "__main__":
    main()
