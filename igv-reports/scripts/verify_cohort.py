#!/usr/bin/env python3
"""verify_cohort.py — cohort-level structural verifier for create_report runs.

Author: Samuel Ahuno
Purpose:
  Catches sample-to-HTML mismatches in cohort mode. Per-sample verification
  (verify_report.py) confirms each HTML is internally consistent, but it
  has no notion of WHICH sample an HTML *should* belong to. This verifier
  re-reads the samplesheet and cross-checks every HTML against the row that
  produced it, plus scans for cross-sample contamination.

Threat model — failure modes this catches that per-sample verify cannot:
  * Wrong BAM embedded under right filename (samplesheet typo, copy-paste).
  * Tumor/normal slot swap.
  * Missing HTML for a samplesheet row (cohort loop silently skipped).
  * Index.html lying — links to a sample that doesn't exist, or omits one.
  * Sample-2's BAM accidentally winding up inside sample-1's HTML.

Checks emitted (per sample, plus two cohort-global rows tagged sample="*"):
  Per-sample (delegated to verify_report.py for the structural ones):
    * html_exists, html_min_size, region_count, region_coords,
      region_sessions, tracks_present  -- run verify_report.py against
      each sample's HTML using that sample's row as input
  Cohort-specific (added here):
    C2 sample_tracks_match           -- the HTML's session contains every
                                        track basename declared in this row
    C3 no_cross_sample_contamination -- the HTML's session contains NO
                                        basename that belongs to another
                                        row's track columns but not this
                                        row (default-track basenames from
                                        databases_config.yaml are excluded)
    C4 sample_id_embedded            -- the `sample` column value appears in
                                        the HTML's <title> or filename
  Cohort-global (one row each, sample='*'):
    C1 cohort_html_coverage          -- every samplesheet sample has exactly
                                        one matching HTML; flag missing+extras
    C5 index_consistency             -- index.html (if present) links exactly
                                        the samplesheet sample set; each link
                                        target exists and is non-empty

Output:
  TSV with columns: sample / check / status / observed / expected / details
  (also printed to stdout). Optional --summary <path>.md emits a one-page
  rollup: total samples, PASS/FAIL counts per check, contamination incidents
  listed by sample.

Exit code: 0, or 1 if --fail-on-fail is set and any row is FAIL.

Typical use (auto-invoked by build_igvreports.py --samplesheet, but can be
run standalone too):

  python verify_cohort.py \\
      --samplesheet samplesheet.tsv \\
      --reports-dir results/<run>/reports/ \\
      --genome hg38 \\
      --out results/<run>/reports/cohort_verify.tsv \\
      --summary results/<run>/reports/cohort_verify.summary.md \\
      --fail-on-fail
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

# Same-dir imports — both verify_report.py and build_igvreports.py live here.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_report as vr
import build_igvreports as bir


_DBCONFIG_ENV = os.environ.get("IGV_REPORTS_DB_CONFIG")
DEFAULT_DBCONFIG = Path(_DBCONFIG_ENV) if _DBCONFIG_ENV else None
DEFAULT_TRACK_COLUMNS = ["bam_tumor", "bam_normal", "vcf", "extra_tracks"]


@dataclass
class CohortCheck:
    sample: str          # "*" for cohort-global checks
    name: str
    status: str          # PASS | FAIL | SKIP
    observed: str = ""
    expected: str = ""
    details: str = ""


# ---------------------------------------------------------------------------
# Samplesheet inspection
# ---------------------------------------------------------------------------

def row_track_paths(row: dict, track_columns: list[str]) -> list[str]:
    """Extract all track paths from a samplesheet row. Honors `extra_tracks`
    being a comma-separated list (per build_igvreports.py convention)."""
    paths: list[str] = []
    for col in track_columns:
        val = row.get(col)
        if not val or not val.strip():
            continue
        if col == "extra_tracks":
            paths.extend(p.strip() for p in val.split(",") if p.strip())
        else:
            paths.append(val.strip())
    return paths


def track_labels_of(paths: list[str]) -> set[str]:
    """Return the names igv-reports auto-assigns to positional --tracks for
    these paths. igv-reports strips ONE final suffix (verified against
    create_report 1.16.2 — see verify_report.expected_track_labels)."""
    return {Path(p).stem for p in paths}


# ---------------------------------------------------------------------------
# Cohort-global checks (C1, C5)
# ---------------------------------------------------------------------------

def check_html_coverage(rows: list[dict], reports_dir: Path, genome: str) -> CohortCheck:
    expected_files = {f"{r['sample']}.{genome}.html" for r in rows}
    actual_files = {p.name for p in reports_dir.glob(f"*.{genome}.html")}
    missing = sorted(expected_files - actual_files)
    extras = sorted(actual_files - expected_files)
    if not missing and not extras:
        return CohortCheck(
            "*", "cohort_html_coverage", "PASS",
            observed=f"{len(actual_files)} HTMLs",
            expected=f"{len(expected_files)} HTMLs",
        )
    details = []
    if missing:
        details.append(f"missing: {', '.join(missing[:5])}" + (" ..." if len(missing) > 5 else ""))
    if extras:
        details.append(f"unexpected: {', '.join(extras[:5])}" + (" ..." if len(extras) > 5 else ""))
    return CohortCheck(
        "*", "cohort_html_coverage", "FAIL",
        observed=f"{len(actual_files)} HTMLs",
        expected=f"{len(expected_files)} HTMLs",
        details="; ".join(details),
    )


def check_index_consistency(rows: list[dict], reports_dir: Path) -> CohortCheck:
    index = reports_dir / "index.html"
    if not index.exists():
        return CohortCheck(
            "*", "index_consistency", "SKIP",
            details=f"no {index.name} present (cohort write_index() not invoked)",
        )
    text = index.read_text()
    # build_igvreports.write_index() emits <li><a href="<file>">SAMPLE</a></li>.
    # Match <a href="..."> ... </a> and pull both the href and the link text.
    found: dict[str, str] = {}  # sample -> href
    for m in re.finditer(r'<a href="([^"]+)">([^<]+)</a>', text):
        href, label = m.group(1), m.group(2).strip()
        found[label] = href

    expected_samples = {r["sample"] for r in rows}
    indexed_samples = set(found.keys())
    missing = sorted(expected_samples - indexed_samples)
    extras = sorted(indexed_samples - expected_samples)
    broken_links = []
    for sample, href in found.items():
        target = reports_dir / href
        if not target.exists() or target.stat().st_size < 1024:
            broken_links.append(f"{sample}->{href}")

    if not missing and not extras and not broken_links:
        return CohortCheck(
            "*", "index_consistency", "PASS",
            observed=f"{len(found)} links",
            expected=f"{len(expected_samples)} samples",
        )
    details = []
    if missing:
        details.append(f"missing from index: {', '.join(missing[:5])}")
    if extras:
        details.append(f"unexpected in index: {', '.join(extras[:5])}")
    if broken_links:
        details.append(f"broken: {', '.join(broken_links[:5])}")
    return CohortCheck(
        "*", "index_consistency", "FAIL",
        observed=f"{len(found)} links",
        expected=f"{len(expected_samples)} samples",
        details="; ".join(details),
    )


# ---------------------------------------------------------------------------
# Per-sample checks (delegate to verify_report + add C2, C3, C4)
# ---------------------------------------------------------------------------

def per_sample_structural(sample: str, html_path: Path, sites_path: Path,
                          tracks: list[str], min_size_mb: float) -> list[CohortCheck]:
    """Run verify_report.py's 6 structural checks against one sample's HTML."""
    out: list[CohortCheck] = []
    out.append(_wrap(sample, vr.check_html_exists(html_path)))
    if not html_path.is_file():
        for n in ("html_min_size", "region_count", "region_coords",
                  "region_sessions", "tracks_present"):
            out.append(CohortCheck(sample, n, "SKIP", details="HTML missing"))
        return out
    out.append(_wrap(sample, vr.check_html_min_size(html_path, min_size_mb)))
    if not sites_path.exists():
        for n in ("region_count", "region_coords", "region_sessions", "tracks_present"):
            out.append(CohortCheck(sample, n, "SKIP", details=f"sites BED missing: {sites_path}"))
        return out
    html_text = html_path.read_text()
    table_json = vr.parse_table_json(html_text)
    session_dict = vr.parse_session_dictionary(html_text)
    bed_rows = vr.load_sites_bed(sites_path)
    out.append(_wrap(sample, vr.check_region_count(bed_rows, table_json)))
    out.append(_wrap(sample, vr.check_region_coords(bed_rows, table_json)))
    out.append(_wrap(sample, vr.check_region_sessions(table_json, session_dict)))
    labels = vr.expected_track_labels(tracks, None)
    out.append(_wrap(sample, vr.check_tracks_present(session_dict, labels)))
    return out


def _wrap(sample: str, c: vr.Check) -> CohortCheck:
    return CohortCheck(sample, c.name, c.status, c.observed, c.expected, c.details)


def session_track_names(html_path: Path) -> set[str]:
    """Decode the first sessionDictionary entry and return its track names.
    Returns an empty set on any decode failure."""
    if not html_path.is_file():
        return set()
    text = html_path.read_text()
    sd = vr.parse_session_dictionary(text)
    if not sd:
        return set()
    sample_key = sorted(sd.keys())[0]
    session = vr.decode_session_entry(sd[sample_key])
    if session is None:
        return set()
    return {t.get("name") for t in session.get("tracks", []) if t.get("name")}


def check_sample_tracks_match(sample: str, html_path: Path, row_tracks: list[str]) -> CohortCheck:
    """C2: each track-stem declared in this sample's row appears as a track
    name in this HTML's session. (igv-reports auto-names positional tracks
    by Path.stem — see verify_report.py's expected_track_labels rationale.)"""
    if not html_path.is_file():
        return CohortCheck(sample, "sample_tracks_match", "SKIP", details="HTML missing")
    expected = sorted(track_labels_of(row_tracks))
    if not expected:
        return CohortCheck(sample, "sample_tracks_match", "SKIP",
                           details="no track paths in samplesheet row")
    names = session_track_names(html_path)
    misses = [b for b in expected if b not in names]
    if misses:
        return CohortCheck(
            sample, "sample_tracks_match", "FAIL",
            observed=f"{len(expected) - len(misses)}/{len(expected)} found",
            expected=f"{len(expected)}/{len(expected)} found",
            details="missing: " + ", ".join(misses[:5]) + (" ..." if len(misses) > 5 else ""),
        )
    return CohortCheck(
        sample, "sample_tracks_match", "PASS",
        observed=f"{len(expected)}/{len(expected)} found",
    )


def check_no_cross_sample_contamination(
    sample: str,
    html_path: Path,
    this_row_labels: set[str],
    other_rows_labels: set[str],
    allow_list: set[str],
) -> CohortCheck:
    """C3: HTML must not contain any track-name label that belongs to OTHER
    samplesheet rows but not this one and not the default-track allow list.
    Labels are Path.stem (igv-reports's auto-naming for positional tracks)."""
    if not html_path.is_file():
        return CohortCheck(sample, "no_cross_sample_contamination", "SKIP", details="HTML missing")
    suspicious = (other_rows_labels - this_row_labels) - allow_list
    if not suspicious:
        return CohortCheck(
            sample, "no_cross_sample_contamination", "PASS",
            observed="0 suspect labels in scope",
        )
    names = session_track_names(html_path)
    incidents = sorted([b for b in suspicious if b in names])
    if not incidents:
        return CohortCheck(
            sample, "no_cross_sample_contamination", "PASS",
            observed=f"{len(suspicious)} other-sample labels scanned, 0 found",
        )
    return CohortCheck(
        sample, "no_cross_sample_contamination", "FAIL",
        observed=f"{len(incidents)} contamination incidents",
        details="found: " + ", ".join(incidents[:5]) + (" ..." if len(incidents) > 5 else ""),
    )


def check_sample_id_embedded(sample: str, html_path: Path) -> CohortCheck:
    """C4: the sample id appears in the HTML's embedded <title>.

    Filename is intentionally NOT checked. The filename is what the cohort
    loop named the file; the title is what `create_report --title` baked
    INTO the HTML at render time. For swap detection, only the title is a
    real signal — a copy-paste of sample_2.html over sample_1.html leaves
    the filename as `sample_1.hg38.html` but the title still says
    `sample_2 (hg38)`. Build_igvreports.py's default title pattern is
    `<sample> (<genome>)`, so this works out of the box.

    If --title is overridden and omits the sample id, this check will FAIL
    — which is the right behavior for a verifier that doesn't know the
    user's intent."""
    if not html_path.is_file():
        return CohortCheck(sample, "sample_id_embedded", "SKIP", details="HTML missing")
    # Read just the head so we don't scan 25 MB for a string.
    head = html_path.read_text()[:16384]
    m = re.search(r"<title>([^<]*)</title>", head, flags=re.IGNORECASE)
    if not m:
        return CohortCheck(
            sample, "sample_id_embedded", "SKIP",
            details="no <title> tag in HTML head; cannot verify",
        )
    title = m.group(1)
    if sample in title:
        return CohortCheck(sample, "sample_id_embedded", "PASS",
                           observed=f"in <title>: {title!r}")
    return CohortCheck(
        sample, "sample_id_embedded", "FAIL",
        observed=f"title={title!r}",
        details=f"sample id {sample!r} not in <title> — likely a swap or wrong --title",
    )


# ---------------------------------------------------------------------------
# Allow-list (default tracks resolved from databases_config.yaml)
# ---------------------------------------------------------------------------

def resolve_default_track_labels(db_config: Path, genome: str) -> set[str]:
    """Reuse the driver's logic so the allow-list stays in sync with what was
    actually loaded. Returns Path.stem of each default track (matches igv-
    reports's auto-naming convention — see track_labels_of)."""
    import logging
    log = logging.getLogger("verify_cohort.allow_list_probe")
    log.addHandler(logging.NullHandler())
    cfg = bir.load_db_config(db_config)
    canon = bir.resolve_genome(genome)
    try:
        paths = bir.resolve_default_tracks(cfg, canon, log)
    except SystemExit:
        # genome not in db_config — fail open with an empty allow-list; the
        # contamination check will then be over-conservative, never under.
        return set()
    return {Path(p).stem for p in paths}


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_tsv(checks: list[CohortCheck], out: Path | None) -> None:
    lines = ["sample\tcheck\tstatus\tobserved\texpected\tdetails"]
    for c in checks:
        lines.append("\t".join((c.sample, c.name, c.status, c.observed, c.expected, c.details)))
    text = "\n".join(lines) + "\n"
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    sys.stdout.write(text)


def write_summary(checks: list[CohortCheck], rows: list[dict], out: Path) -> None:
    by_status = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    by_check: dict[str, dict[str, int]] = {}
    fail_rows = []
    for c in checks:
        by_status[c.status] = by_status.get(c.status, 0) + 1
        by_check.setdefault(c.name, {"PASS": 0, "FAIL": 0, "SKIP": 0})[c.status] += 1
        if c.status == "FAIL":
            fail_rows.append(c)

    n_samples = len(rows)
    lines = []
    lines.append(f"# Cohort verification summary\n")
    lines.append(f"- samples: **{n_samples}**")
    lines.append(f"- total checks: {sum(by_status.values())} (PASS={by_status['PASS']}, FAIL={by_status['FAIL']}, SKIP={by_status['SKIP']})")
    lines.append("")
    lines.append("## Per-check totals")
    lines.append("")
    lines.append("| check | PASS | FAIL | SKIP |")
    lines.append("|---|---:|---:|---:|")
    for check_name in sorted(by_check):
        s = by_check[check_name]
        lines.append(f"| {check_name} | {s['PASS']} | {s['FAIL']} | {s['SKIP']} |")
    lines.append("")
    if fail_rows:
        lines.append("## Failures")
        lines.append("")
        lines.append("| sample | check | observed | expected | details |")
        lines.append("|---|---|---|---|---|")
        for c in fail_rows:
            lines.append(f"| {c.sample} | {c.name} | {c.observed} | {c.expected} | {c.details} |")
    else:
        lines.append("## Failures\n\nNone — cohort verified clean.\n")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True, help="TSV that drove the cohort build (same one passed to build_igvreports.py --samplesheet)")
    ap.add_argument("--reports-dir", required=True, help="dir containing <sample>.<genome>.html files (and optional index.html)")
    ap.add_argument("--genome", required=True, help="genome tag (hg38 | mm10 | mm39 | t2t | GRCh37)")
    ap.add_argument("--db-config", default=(str(DEFAULT_DBCONFIG) if DEFAULT_DBCONFIG else None), help="databases_config.yaml to resolve default-track allow-list (optional)")
    ap.add_argument(
        "--track-columns", nargs="*", default=DEFAULT_TRACK_COLUMNS,
        help=f"samplesheet columns containing track paths (default: {DEFAULT_TRACK_COLUMNS}). "
             "`extra_tracks` is parsed comma-separated if present.",
    )
    ap.add_argument("--min-size-mb", type=float, default=0.5, help="per-sample HTML min size (passed through to verify_report)")
    ap.add_argument("--out", help="write the TSV report here in addition to stdout")
    ap.add_argument("--summary", help="write a one-page markdown rollup here")
    ap.add_argument("--fail-on-fail", action="store_true", help="exit nonzero if any check is FAIL")
    args = ap.parse_args()

    samplesheet = Path(args.samplesheet)
    reports_dir = Path(args.reports_dir)
    if not samplesheet.exists():
        raise SystemExit(f"ERROR: samplesheet not found: {samplesheet}")
    if not reports_dir.is_dir():
        raise SystemExit(f"ERROR: reports-dir not found: {reports_dir}")

    rows = bir.parse_samplesheet(samplesheet)
    if not rows:
        raise SystemExit(f"ERROR: samplesheet has no data rows: {samplesheet}")

    allow_list = (
        resolve_default_track_labels(Path(args.db_config), args.genome)
        if args.db_config else set()
    )

    # Pre-compute track-label sets per sample for the contamination check.
    # Labels are Path.stem of each track path, matching igv-reports's auto-
    # naming (see track_labels_of).
    per_sample_labels: dict[str, set[str]] = {
        r["sample"]: track_labels_of(row_track_paths(r, args.track_columns)) for r in rows
    }
    all_labels = set().union(*per_sample_labels.values()) if per_sample_labels else set()

    checks: list[CohortCheck] = []
    # C1 cohort_html_coverage
    checks.append(check_html_coverage(rows, reports_dir, args.genome))

    # Per-sample: 6 structural (verify_report) + C2 + C3 + C4
    for r in rows:
        sample = r["sample"]
        html_path = reports_dir / f"{sample}.{args.genome}.html"
        sites_path = Path(r["sites_bed"])
        tracks = row_track_paths(r, args.track_columns)

        checks.extend(per_sample_structural(sample, html_path, sites_path, tracks, args.min_size_mb))
        checks.append(check_sample_tracks_match(sample, html_path, tracks))

        this_labels = per_sample_labels[sample]
        other_labels = all_labels - this_labels
        checks.append(check_no_cross_sample_contamination(sample, html_path, this_labels, other_labels, allow_list))
        checks.append(check_sample_id_embedded(sample, html_path))

    # C5 index_consistency
    checks.append(check_index_consistency(rows, reports_dir))

    out_path = Path(args.out) if args.out else None
    write_tsv(checks, out_path)
    if args.summary:
        write_summary(checks, rows, Path(args.summary))

    if args.fail_on_fail and any(c.status == "FAIL" for c in checks):
        sys.exit(1)


if __name__ == "__main__":
    main()
