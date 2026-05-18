#!/usr/bin/env python3
"""verify_anchors.py — content verifier for create_report HTMLs.

Author: Samuel Ahuno
Purpose:
  The structural verifier (verify_report / verify_cohort) confirms the HTML
  *says* the right thing: region count, coords, track names. It cannot
  confirm the embedded BAM slices actually contain the data they claim to.
  Failure modes it misses:

    1. Sample swap — track name says `p17424_1.sorted` but the slice was
       cut from `p17424_3.sorted.bam` (cohort loop wired the wrong path;
       Path.stem matched and the structural check passed).
    2. Silent empty slice — region rendered, slice is 0 reads (failed
       index, BAM corruption, coords outside coverage).
    3. Regression across create_report versions — flanking/slicing logic
       changes silently between releases.

  This verifier closes the gap by re-running `samtools view -c` against
  both the source BAM (at generate time) and the embedded slice (at
  verify time), then comparing counts.

  Anchor TSV format (`#`-prefixed header, lab BED-output convention):

    #sample	track_name	chrom	start	end	expected	tolerance	min	max	notes

  - `tolerance` and `min`/`max` are mutually exclusive per row; if `min`
    or `max` is non-empty it wins. Blank tolerance falls back to
    --tolerance flag default (0.05).
  - `expected` is the count from `samtools view -c -F 1536 source.bam
    chrom:start-end` at generate time. Generate writes it; verify reads it.

Subcommands:
  generate      — walk (sample × region) grid, count reads from source BAMs,
                  write an anchors.tsv that becomes a regression fixture.
  verify        — given one HTML + anchors.tsv, decode each anchor's BAM
                  slice and count it, compare to expected.
  verify-cohort — apply `verify` across all HTMLs in a cohort.

Container resolution (samtools):
  1. --samtools-sif PATH
  2. $SAMTOOLS_SIF_DEFAULT env var
  3. `samtools` on PATH
  4. Hard error

Typical use:
  # at build time, freeze the regression fixture
  python verify_anchors.py generate \\
      --samplesheet sheet.tsv \\
      --sites sites.hg38.bed \\
      --out anchors.hg38.tsv

  # any time after, audit a built HTML
  python verify_anchors.py verify \\
      --html report.hg38.html \\
      --anchors anchors.hg38.tsv \\
      --out verify_anchors.tsv \\
      --fail-on-fail

  # cohort-wide
  python verify_anchors.py verify-cohort \\
      --samplesheet sheet.tsv \\
      --reports-dir results/<run>/reports/ \\
      --genome hg38 \\
      --anchors anchors.hg38.tsv \\
      --out cohort_verify_anchors.tsv
"""

from __future__ import annotations

import argparse
import base64
import dataclasses
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# Same-dir imports — reuse verify_report's HTML parser helpers.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_report as vr


_SAMTOOLS_SIF_ENV = os.environ.get("SAMTOOLS_SIF_DEFAULT")
DEFAULT_SAMTOOLS_SIF = Path(_SAMTOOLS_SIF_ENV) if _SAMTOOLS_SIF_ENV else None


def _apptainer_bind_args() -> list[str]:
    """Conditional `--bind <path>` tokens, matching build_igvreports.py.

    Source: $IGV_REPORTS_BIND (colon-separated). Empty / unset = no binds.
    Paths that don't exist are silently skipped."""
    raw = os.environ.get("IGV_REPORTS_BIND", "")
    if not raw:
        return []
    tokens: list[str] = []
    for p in raw.split(":"):
        if p and Path(p).exists():
            tokens.extend(["--bind", p])
    return tokens
# Match igv-reports BamReader default exclude flag (rules out PCR/optical
# duplicates and supplementary alignments — see igv_reports/bam.py).
EXCLUDE_FLAGS = "1536"
DEFAULT_TOLERANCE = 0.05
ANCHOR_HEADER = [
    "sample", "track_name", "chrom", "start", "end",
    "expected", "tolerance", "min", "max", "notes",
]


@dataclasses.dataclass
class AnchorRow:
    sample: str
    track_name: str
    chrom: str
    start: int
    end: int
    expected: int
    tolerance: str = ""   # blank => fall back to --tolerance flag
    min_count: str = ""   # blank => not used
    max_count: str = ""   # blank => not used
    notes: str = ""

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclasses.dataclass
class AnchorCheck:
    sample: str
    track_name: str
    region: str
    status: str           # PASS | FAIL | SKIP
    observed: str = ""
    expected: str = ""
    details: str = ""


# ---------------------------------------------------------------------------
# samtools resolution
# ---------------------------------------------------------------------------

def resolve_samtools(sif: Path | None) -> list[str]:
    """Return a samtools command prefix (list of argv tokens).

    Priority: --samtools-sif → $SAMTOOLS_SIF → $SAMTOOLS_SIF_DEFAULT → PATH.
    Falling back to PATH is fine for non-HPC use; on HPC, prefer a SIF to
    avoid the NFS conda cold-start tax.
    """
    candidate = sif
    if candidate is None:
        env = os.environ.get("SAMTOOLS_SIF")
        if env:
            candidate = Path(env)
    if candidate is None and DEFAULT_SAMTOOLS_SIF is not None and DEFAULT_SAMTOOLS_SIF.exists():
        candidate = DEFAULT_SAMTOOLS_SIF
    if candidate is not None:
        if not candidate.exists():
            raise SystemExit(f"ERROR: samtools SIF not found: {candidate}")
        return [
            "singularity", "exec", "--cleanenv", *_apptainer_bind_args(),
            str(candidate), "samtools",
        ]
    path_sam = shutil.which("samtools")
    if path_sam:
        sys.stderr.write(
            f"[verify_anchors] WARNING: falling back to PATH samtools at {path_sam}; "
            "SIF preferred for HPC cold-start cost (rules/apptainer_vs_conda.md)\n"
        )
        return [path_sam]
    raise SystemExit(
        "ERROR: no samtools found. Provide --samtools-sif, set $SAMTOOLS_SIF, "
        "or install samtools on PATH."
    )


def samtools_count(samtools_cmd: list[str], bam: Path, region: str) -> int:
    """Run `samtools view -c -F 1536 <bam> <region>` and return the count."""
    proc = subprocess.run(
        samtools_cmd + ["view", "-c", "-F", EXCLUDE_FLAGS, str(bam), region],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"samtools view -c failed (exit {proc.returncode}) for {bam} {region}: "
            f"{proc.stderr.strip()}"
        )
    return int(proc.stdout.strip())


def samtools_index(samtools_cmd: list[str], bam: Path) -> None:
    """Run `samtools index <bam>`."""
    proc = subprocess.run(
        samtools_cmd + ["index", str(bam)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"samtools index failed (exit {proc.returncode}) for {bam}: "
            f"{proc.stderr.strip()}"
        )


# ---------------------------------------------------------------------------
# anchors.tsv I/O
# ---------------------------------------------------------------------------

def write_anchors(anchors: list[AnchorRow], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#" + "\t".join(ANCHOR_HEADER)]
    for a in anchors:
        lines.append("\t".join([
            a.sample, a.track_name, a.chrom, str(a.start), str(a.end),
            str(a.expected), a.tolerance, a.min_count, a.max_count, a.notes,
        ]))
    out.write_text("\n".join(lines) + "\n")


def load_anchors(path: Path) -> list[AnchorRow]:
    if not path.exists():
        raise SystemExit(f"ERROR: anchors TSV not found: {path}")
    rows: list[AnchorRow] = []
    with path.open() as fh:
        header: list[str] | None = None
        for i, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").split("\t")
                continue
            if header is None:
                raise SystemExit(f"{path}:{i}: data row before header — anchors TSV needs a `#`-prefixed header")
            cols = line.split("\t")
            if len(cols) < len(header):
                cols += [""] * (len(header) - len(cols))
            d = dict(zip(header, cols))
            try:
                # Validate numeric optional fields at load time so a mis-tabbed
                # row fails here, not deep inside decide_status() with a
                # confusing 'could not convert' on the notes value.
                tolerance = (d.get("tolerance", "") or "").strip()
                if tolerance:
                    float(tolerance)
                min_count = (d.get("min", "") or "").strip()
                if min_count:
                    int(min_count)
                max_count = (d.get("max", "") or "").strip()
                if max_count:
                    int(max_count)
                rows.append(AnchorRow(
                    sample=d["sample"],
                    track_name=d["track_name"],
                    chrom=d["chrom"],
                    start=int(d["start"]),
                    end=int(d["end"]),
                    expected=int(d["expected"]),
                    tolerance=tolerance,
                    min_count=min_count,
                    max_count=max_count,
                    notes=d.get("notes", "") or "",
                ))
            except (KeyError, ValueError) as e:
                raise SystemExit(
                    f"{path}:{i}: malformed anchor row: {e}\n"
                    f"  row was: {cols!r}\n"
                    f"  expected columns: {ANCHOR_HEADER}\n"
                    f"  hint: TSV reader requires explicit tab separation — "
                    "if you generate the row with awk, pass `-F'\\t'`."
                )
    return rows


# ---------------------------------------------------------------------------
# samplesheet → (sample, track_path) iteration (shared with build_igvreports)
# ---------------------------------------------------------------------------

def parse_samplesheet(path: Path) -> list[dict]:
    """Mirror build_igvreports.parse_samplesheet without importing it (avoids
    pulling in PyYAML for code paths that don't need it)."""
    rows: list[dict] = []
    with path.open() as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for ln in fh:
            cols = ln.rstrip("\n").split("\t")
            if not cols or not cols[0].strip():
                continue
            rows.append(dict(zip(header, cols)))
    return rows


def sample_bam_paths(row: dict) -> list[tuple[str, Path]]:
    """Return [(track_name, bam_path), ...] for the BAM columns in a row.
    track_name = Path.stem (matches igv-reports' positional auto-naming —
    see verify_report.expected_track_labels)."""
    out: list[tuple[str, Path]] = []
    for col in ("bam_tumor", "bam_normal"):
        v = row.get(col)
        if v and v.strip():
            p = Path(v.strip())
            out.append((p.stem, p))
    extras = row.get("extra_tracks") or ""
    for entry in extras.split(","):
        entry = entry.strip()
        if entry.endswith(".bam") or entry.endswith(".cram"):
            p = Path(entry)
            out.append((p.stem, p))
    return out


# ---------------------------------------------------------------------------
# Slice extraction from embedded session
# ---------------------------------------------------------------------------

_DATA_URL_RE = re.compile(r"data:[^;]+;base64,(.+)", flags=re.DOTALL)


def decode_track_slice(track_url: str, dest: Path) -> Path:
    """Decode a track's `data:...;base64,...` URL, write bytes to `dest`.

    Per igv_reports/datauri.py: BAM slices come back from pysam.view as
    bytes starting with BGZF magic (0x1f 0x8b), so igv-reports tags them
    as `data:application/gzip;base64,...`. We accept any data: URL with a
    base64 payload — the bytes are what matters, not the declared mediatype.
    """
    m = _DATA_URL_RE.match(track_url)
    if not m:
        raise ValueError("track url is not a data: base64 URL")
    raw = base64.b64decode(m.group(1))
    dest.write_bytes(raw)
    return dest


def locate_session_entry(
    session_dict: dict, table_json: dict, chrom: str, start: int, end: int,
) -> tuple[str, dict | None, str]:
    """Locate the session entry for an anchor's (chrom, start+1, end).

    Returns (outcome, session_or_none, detail) where outcome is one of:
      'absent'   — no tableJson row matches this region → caller should SKIP
                   (anchor lists a region the HTML never rendered)
      'broken'   — row matched but session missing/undecodable → caller FAILs
                   (structural inconsistency or HTML corruption)
      'ok'       — session decoded; second element is the dict
    HTML stores 1-based start (per verify_report comment); BED is 0-based.
    """
    headers = table_json.get("headers", [])
    try:
        col_chrom = headers.index("Chrom")
        col_start = headers.index("Start")
        col_end = headers.index("End")
    except ValueError as e:
        return ("broken", None, f"tableJson missing expected column: {e}")
    rows = table_json.get("rows", [])
    want = (chrom, start + 1, end)
    for idx, row in enumerate(rows):
        if (row[col_chrom], int(row[col_start]), int(row[col_end])) == want:
            data_url = session_dict.get(str(idx))
            if data_url is None:
                return ("broken", None, f"sessionDictionary has no entry for row index {idx}")
            session = vr.decode_session_entry(data_url)
            if session is None:
                return ("broken", None, f"session entry {idx} failed to gunzip/decode")
            return ("ok", session, "")
    return ("absent", None, f"no tableJson row matched ({chrom}, {start+1}, {end})")


def find_track(session: dict, track_name: str) -> dict | None:
    for t in session.get("tracks", []):
        if t.get("name") == track_name:
            return t
    return None


# ---------------------------------------------------------------------------
# Status decision
# ---------------------------------------------------------------------------

def decide_status(anchor: AnchorRow, observed: int, default_tol: float) -> tuple[str, str]:
    """Return (status, details). min/max wins over tolerance when present."""
    if anchor.min_count or anchor.max_count:
        bounds_ok = True
        bits = []
        if anchor.min_count:
            ok = observed >= int(anchor.min_count)
            bits.append(f"min={anchor.min_count} {'OK' if ok else 'FAIL'}")
            bounds_ok = bounds_ok and ok
        if anchor.max_count:
            ok = observed <= int(anchor.max_count)
            bits.append(f"max={anchor.max_count} {'OK' if ok else 'FAIL'}")
            bounds_ok = bounds_ok and ok
        return ("PASS" if bounds_ok else "FAIL"), "; ".join(bits)
    tol = float(anchor.tolerance) if anchor.tolerance else default_tol
    if anchor.expected == 0:
        ok = observed == 0
        return ("PASS" if ok else "FAIL"), f"expected=0, observed={observed}"
    diff_ratio = abs(observed - anchor.expected) / anchor.expected
    ok = diff_ratio <= tol
    return ("PASS" if ok else "FAIL"), f"diff_ratio={diff_ratio:.3f} (tol={tol:.3f})"


# ---------------------------------------------------------------------------
# Subcommand: generate
# ---------------------------------------------------------------------------

def cmd_generate(args: argparse.Namespace) -> None:
    samtools_cmd = resolve_samtools(Path(args.samtools_sif) if args.samtools_sif else None)
    rows = parse_samplesheet(Path(args.samplesheet))
    bed_rows = vr.load_sites_bed(Path(args.sites))
    if not rows:
        raise SystemExit("ERROR: samplesheet has no data rows")
    if not bed_rows:
        raise SystemExit("ERROR: sites BED has no data rows")

    anchors: list[AnchorRow] = []
    for row in rows:
        sample = row["sample"]
        bams = sample_bam_paths(row)
        if not bams:
            sys.stderr.write(f"[generate] {sample}: no BAM tracks in row — skipping\n")
            continue
        for track_name, bam in bams:
            if not bam.exists():
                sys.stderr.write(f"[generate] {sample}/{track_name}: BAM missing: {bam}\n")
                continue
            for b in bed_rows:
                region = f"{b['chrom']}:{b['start']}-{b['end']}"
                try:
                    count = samtools_count(samtools_cmd, bam, region)
                except RuntimeError as e:
                    sys.stderr.write(f"[generate] {sample}/{track_name} {region}: {e}\n")
                    continue
                anchors.append(AnchorRow(
                    sample=sample, track_name=track_name,
                    chrom=b["chrom"], start=b["start"], end=b["end"],
                    expected=count, notes=b["name"] or "",
                ))
                sys.stderr.write(f"[generate] {sample}/{track_name} {region}: {count}\n")

    out = Path(args.out)
    write_anchors(anchors, out)
    sys.stderr.write(f"[generate] wrote {len(anchors)} anchors -> {out}\n")


# ---------------------------------------------------------------------------
# Subcommand: verify (single HTML)
# ---------------------------------------------------------------------------

def verify_one_html(
    html_path: Path, anchors: list[AnchorRow], samtools_cmd: list[str],
    default_tol: float,
) -> list[AnchorCheck]:
    """Verify all anchors against one HTML. Anchors whose track_name doesn't
    appear in the HTML are SKIPped (cohort verify-cohort filters by sample,
    so this function trusts the caller passed the right anchor subset)."""
    checks: list[AnchorCheck] = []
    if not html_path.is_file():
        for a in anchors:
            checks.append(AnchorCheck(
                a.sample, a.track_name, a.region, "SKIP",
                details=f"HTML missing: {html_path}",
            ))
        return checks
    html_text = html_path.read_text()
    table_json = vr.parse_table_json(html_text)
    session_dict = vr.parse_session_dictionary(html_text)
    if table_json is None or session_dict is None:
        for a in anchors:
            checks.append(AnchorCheck(
                a.sample, a.track_name, a.region, "FAIL",
                details="tableJson or sessionDictionary missing from HTML",
            ))
        return checks
    with tempfile.TemporaryDirectory(prefix="verify_anchors_") as td:
        tmp = Path(td)
        for a in anchors:
            outcome, session, locate_detail = locate_session_entry(
                session_dict, table_json, a.chrom, a.start, a.end,
            )
            if outcome == "absent":
                checks.append(AnchorCheck(
                    a.sample, a.track_name, a.region, "SKIP",
                    details=locate_detail,
                ))
                continue
            if outcome == "broken":
                checks.append(AnchorCheck(
                    a.sample, a.track_name, a.region, "FAIL",
                    expected=str(a.expected),
                    details=locate_detail,
                ))
                continue
            assert session is not None  # outcome == "ok"
            track = find_track(session, a.track_name)
            if track is None:
                checks.append(AnchorCheck(
                    a.sample, a.track_name, a.region, "SKIP",
                    details=f"track '{a.track_name}' not in HTML session",
                ))
                continue
            url = track.get("url", "")
            slice_path = tmp / f"{a.sample}__{a.track_name}__{a.chrom}_{a.start}_{a.end}.bam"
            try:
                decode_track_slice(url, slice_path)
                samtools_index(samtools_cmd, slice_path)
                observed = samtools_count(samtools_cmd, slice_path, a.region)
            except (ValueError, RuntimeError) as e:
                checks.append(AnchorCheck(
                    a.sample, a.track_name, a.region, "FAIL",
                    expected=str(a.expected),
                    details=f"slice decode/count failed: {e}",
                ))
                continue
            status, details = decide_status(a, observed, default_tol)
            checks.append(AnchorCheck(
                a.sample, a.track_name, a.region, status,
                observed=str(observed), expected=str(a.expected),
                details=details,
            ))
    return checks


def cmd_verify(args: argparse.Namespace) -> None:
    samtools_cmd = resolve_samtools(Path(args.samtools_sif) if args.samtools_sif else None)
    anchors = load_anchors(Path(args.anchors))
    checks = verify_one_html(Path(args.html), anchors, samtools_cmd, args.tolerance)
    write_checks(checks, Path(args.out) if args.out else None)
    if args.fail_on_fail and any(c.status == "FAIL" for c in checks):
        sys.exit(1)


# ---------------------------------------------------------------------------
# Subcommand: verify-cohort
# ---------------------------------------------------------------------------

def cmd_verify_cohort(args: argparse.Namespace) -> None:
    samtools_cmd = resolve_samtools(Path(args.samtools_sif) if args.samtools_sif else None)
    anchors = load_anchors(Path(args.anchors))
    rows = parse_samplesheet(Path(args.samplesheet))
    reports_dir = Path(args.reports_dir)
    genome = args.genome

    # Group anchors by sample for per-HTML filtering.
    by_sample: dict[str, list[AnchorRow]] = {}
    for a in anchors:
        by_sample.setdefault(a.sample, []).append(a)

    all_checks: list[AnchorCheck] = []
    for row in rows:
        sample = row["sample"]
        html_path = reports_dir / f"{sample}.{genome}.html"
        sample_anchors = by_sample.get(sample, [])
        if not sample_anchors:
            all_checks.append(AnchorCheck(
                sample, "*", "*", "SKIP",
                details="no anchors for this sample in anchors.tsv",
            ))
            continue
        all_checks.extend(verify_one_html(html_path, sample_anchors, samtools_cmd, args.tolerance))

    # Surface anchor samples that don't match any samplesheet row.
    samplesheet_samples = {r["sample"] for r in rows}
    anchor_orphans = sorted(set(by_sample.keys()) - samplesheet_samples)
    for s in anchor_orphans:
        all_checks.append(AnchorCheck(
            s, "*", "*", "SKIP",
            details="anchor sample not present in samplesheet",
        ))

    write_checks(all_checks, Path(args.out) if args.out else None)
    if args.fail_on_fail and any(c.status == "FAIL" for c in all_checks):
        sys.exit(1)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_checks(checks: list[AnchorCheck], out: Path | None) -> None:
    lines = ["sample\ttrack_name\tregion\tstatus\tobserved\texpected\tdetails"]
    for c in checks:
        lines.append("\t".join([
            c.sample, c.track_name, c.region, c.status,
            c.observed, c.expected, c.details,
        ]))
    text = "\n".join(lines) + "\n"
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    sys.stdout.write(text)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    # generate
    g = sub.add_parser("generate", help="freeze samtools view -c counts into anchors.tsv")
    g.add_argument("--samplesheet", required=True)
    g.add_argument("--sites", required=True)
    g.add_argument("--out", required=True, help="path to write anchors TSV")
    g.add_argument("--samtools-sif", help="explicit samtools SIF path")
    g.set_defaults(func=cmd_generate)

    # verify
    v = sub.add_parser("verify", help="audit one HTML against anchors.tsv")
    v.add_argument("--html", required=True)
    v.add_argument("--anchors", required=True)
    v.add_argument("--out", help="write checks TSV here in addition to stdout")
    v.add_argument("--samtools-sif")
    v.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE,
                   help=f"default ratio tolerance when row tolerance/min/max blank (default {DEFAULT_TOLERANCE})")
    v.add_argument("--fail-on-fail", action="store_true")
    v.set_defaults(func=cmd_verify)

    # verify-cohort
    vc = sub.add_parser("verify-cohort", help="audit all HTMLs in a cohort against anchors.tsv")
    vc.add_argument("--samplesheet", required=True)
    vc.add_argument("--reports-dir", required=True)
    vc.add_argument("--genome", required=True)
    vc.add_argument("--anchors", required=True)
    vc.add_argument("--out")
    vc.add_argument("--samtools-sif")
    vc.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE)
    vc.add_argument("--fail-on-fail", action="store_true")
    vc.set_defaults(func=cmd_verify_cohort)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
