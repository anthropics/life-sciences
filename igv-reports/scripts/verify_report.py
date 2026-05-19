#!/usr/bin/env python3
"""verify_report.py — post-render structural verifier for create_report HTMLs.

Author: Samuel Ahuno
Purpose:
  Validates that a self-contained create_report HTML actually contains what
  its inputs declared. Catches the "silent garbage" failure mode where the
  HTML builds (exit 0, plausible file size) but the content doesn't match the
  user's intent: wrong region count, wrong coordinates, missing tracks, or a
  catastrophic empty render.

Dual role:
  - CLI: `python verify_report.py --html ... --sites ... [--track-config ...]`
  - Library: importable helpers (parse_table_json, parse_session_dictionary,
    decode_session_entry, load_sites_bed, expected_track_labels, the
    `check_*` functions, and the Check dataclass). verify_cohort.py imports
    these to do per-sample checks + add cross-sample assertions.

Checks emitted (one TSV row per check, ordered):
  1. html_exists           Output file is a regular file.
  2. html_min_size         Output >= --min-size-mb (default 0.5 MB).
  3. region_count          tableJson rows count == sites BED data-row count.
  4. region_coords         Each BED row finds a matching (chrom, start+1, end[, name])
                           in the embedded tableJson. BED is 0-based half-open;
                           create_report stores 1-based start in the table.
  5. region_sessions       sessionDictionary has an entry for each tableJson row.
  6. tracks_present        For --track-config <json>: each track's `name` field
                           appears in the decoded session's tracks[].name list.
                           For --tracks <path...>: each path's Path.stem appears
                           in the decoded session's tracks[].name list. igv-
                           reports strips ONE final suffix when auto-naming
                           positional tracks (e.g. `x.5mC.bedgraph` -> `x.5mC`,
                           `gencode.v47.annotation.gff3.gz` -> `gencode.v47.
                           annotation.gff3`). Skipped if neither flag is given.
                           NOTE: --standalone embeds slices as data: URLs, so
                           original URL paths are absent from the session — we
                           match on track NAMES, which are preserved.

Output:
  TSV with columns: check / status / observed / expected / details
  status is one of PASS / FAIL / SKIP.

Exit code:
  0 always, unless --fail-on-fail is set and at least one row is FAIL.

Typical use:
  python verify_report.py \\
      --html report.hg38.html \\
      --sites sites.hg38.bed \\
      --track-config tracks.json \\
      --out verify.tsv \\
      --min-size-mb 1.0 \\
      --fail-on-fail
"""

from __future__ import annotations

import argparse
import base64
import dataclasses
import gzip
import json
import re
import sys
from pathlib import Path


@dataclasses.dataclass
class Check:
    name: str
    status: str  # PASS | FAIL | SKIP
    observed: str = ""
    expected: str = ""
    details: str = ""


# ---------------------------------------------------------------------------
# Sites-BED loader (mirrors create_report's #-skip behavior)
# ---------------------------------------------------------------------------

def load_sites_bed(path: Path) -> list[dict]:
    """Return a list of {chrom, start, end, name} dicts; skips '#' and 'track '."""
    rows: list[dict] = []
    with path.open() as fh:
        for i, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track "):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                raise SystemExit(f"{path}:{i}: BED row has <3 columns")
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError as e:
                raise SystemExit(f"{path}:{i}: non-numeric start/end: {e}")
            rows.append({
                "chrom": cols[0],
                "start": start,
                "end": end,
                "name": cols[3] if len(cols) >= 4 else None,
            })
    return rows


# ---------------------------------------------------------------------------
# HTML extractors
# ---------------------------------------------------------------------------

def _extract_balanced_blob(text: str, anchor: str, opener: str = "{") -> str | None:
    """Find `anchor` in `text`, then return the substring starting at the next
    `opener` and ending at the matched closer. Skips characters inside double-
    quoted strings (with backslash escapes). Returns None if not found."""
    closer = "}" if opener == "{" else "]"
    i = text.find(anchor)
    if i < 0:
        return None
    start = text.find(opener, i)
    if start < 0:
        return None
    depth = 0
    in_str = False
    escape = False
    for j in range(start, len(text)):
        c = text[j]
        if escape:
            escape = False
            continue
        if c == "\\":
            escape = True
            continue
        if c == '"':
            in_str = not in_str
            continue
        if in_str:
            continue
        if c == opener:
            depth += 1
        elif c == closer:
            depth -= 1
            if depth == 0:
                return text[start:j + 1]
    return None


def parse_table_json(html: str) -> dict | None:
    blob = _extract_balanced_blob(html, "tableJson = ", "{")
    if not blob:
        return None
    return json.loads(blob)


def parse_session_dictionary(html: str) -> dict | None:
    blob = _extract_balanced_blob(html, "sessionDictionary = ", "{")
    if not blob:
        return None
    return json.loads(blob)


def decode_session_entry(data_url: str) -> dict | None:
    """A sessionDictionary value looks like 'data:application/gzip;base64,XXXX'.
    Strip the prefix, base64-decode, gunzip, parse JSON. Return the IGV.js
    session dict (or None on any error — failures here are non-fatal)."""
    try:
        m = re.match(r"data:application/gzip;base64,(.+)", data_url, flags=re.DOTALL)
        if not m:
            return None
        raw = base64.b64decode(m.group(1))
        return json.loads(gzip.decompress(raw))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Track-input parser
# ---------------------------------------------------------------------------

def expected_track_labels(tracks: list[str] | None, track_config: Path | None) -> list[str]:
    """Return the track NAMES we expect to see in the embedded igv.js session.

    `--standalone` replaces every track URL with an inlined `data:...` URL after
    slicing, so URL paths are unrecoverable from the embedded session — we have
    to match on track names instead, which the standalone build preserves.

    - For --track-config <json>: use the `name` field of each entry verbatim.
    - For positional --tracks <path...>: use Path(p).stem (igv-reports strips
      ONE final suffix when auto-naming positional tracks — verified 2026-05-16
      against create_report 1.16.2: `colo829bl_PAU59807.5mC.bedgraph` ->
      `colo829bl_PAU59807.5mC`, `gencode.v47.annotation.gff3.gz` ->
      `gencode.v47.annotation.gff3`, `x.bam` -> `x`).
    Empty list means 'check skipped'.
    """
    out: list[str] = []
    if track_config and track_config.exists():
        with track_config.open() as fh:
            cfg = json.load(fh)
        for entry in cfg:
            name = entry.get("name")
            if name:
                out.append(name)
        return out
    if tracks:
        for t in tracks:
            out.append(Path(t).stem)
    return out


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------

def check_html_exists(html: Path) -> Check:
    if html.is_file():
        return Check("html_exists", "PASS", observed=str(html))
    return Check("html_exists", "FAIL", observed=str(html), details="not a regular file")


def check_html_min_size(html: Path, floor_mb: float) -> Check:
    size_mb = html.stat().st_size / 1024 / 1024
    status = "PASS" if size_mb >= floor_mb else "FAIL"
    return Check(
        "html_min_size",
        status,
        observed=f"{size_mb:.2f} MB",
        expected=f">= {floor_mb:.2f} MB",
    )


def check_region_count(bed_rows: list[dict], table_json: dict | None) -> Check:
    if table_json is None:
        return Check("region_count", "FAIL", details="tableJson not found in HTML")
    n_html = len(table_json.get("rows", []))
    n_bed = len(bed_rows)
    return Check(
        "region_count",
        "PASS" if n_html == n_bed else "FAIL",
        observed=str(n_html),
        expected=str(n_bed),
    )


def check_region_coords(bed_rows: list[dict], table_json: dict | None) -> Check:
    """For each BED row, find a matching row in the HTML by (chrom, start+1, end[, name]).
    The HTML stores 1-based start, BED is 0-based half-open."""
    if table_json is None:
        return Check("region_coords", "FAIL", details="tableJson not found")
    headers = table_json.get("headers", [])
    rows = table_json.get("rows", [])
    try:
        col_chrom = headers.index("Chrom")
        col_start = headers.index("Start")
        col_end = headers.index("End")
        col_name = headers.index("Name") if "Name" in headers else None
    except ValueError as e:
        return Check("region_coords", "FAIL", details=f"missing column in tableJson headers: {e}")

    html_set = {
        (r[col_chrom], int(r[col_start]), int(r[col_end])): (r[col_name] if col_name is not None else None)
        for r in rows
    }
    misses: list[str] = []
    for b in bed_rows:
        key = (b["chrom"], b["start"] + 1, b["end"])
        if key not in html_set:
            misses.append(f"{b['chrom']}:{b['start']}-{b['end']}")
            continue
        # If both have a name, names must match.
        if col_name is not None and b["name"] is not None and html_set[key] != b["name"]:
            misses.append(f"{b['chrom']}:{b['start']}-{b['end']} name mismatch (BED={b['name']!r}, HTML={html_set[key]!r})")
    if misses:
        return Check(
            "region_coords", "FAIL",
            observed=f"{len(bed_rows) - len(misses)}/{len(bed_rows)} matched",
            expected=f"{len(bed_rows)}/{len(bed_rows)} matched",
            details="; ".join(misses[:5]) + (" ..." if len(misses) > 5 else ""),
        )
    return Check("region_coords", "PASS", observed=f"{len(bed_rows)}/{len(bed_rows)} matched")


def check_region_sessions(table_json: dict | None, session_dict: dict | None) -> Check:
    if table_json is None or session_dict is None:
        return Check("region_sessions", "FAIL", details="tableJson or sessionDictionary missing")
    n_rows = len(table_json.get("rows", []))
    n_sess = len(session_dict)
    # Sessions are keyed by stringified row index 0..N-1.
    expected_keys = {str(i) for i in range(n_rows)}
    actual_keys = set(session_dict.keys())
    if expected_keys.issubset(actual_keys):
        return Check(
            "region_sessions", "PASS",
            observed=str(n_sess),
            expected=f">={n_rows} (one per row)",
        )
    return Check(
        "region_sessions", "FAIL",
        observed=f"keys={sorted(actual_keys)[:5]}...",
        expected=f"keys 0..{n_rows-1}",
        details=f"missing keys: {sorted(expected_keys - actual_keys)[:5]}",
    )


def check_tracks_present(
    session_dict: dict | None,
    expected_labels: list[str],
) -> Check:
    if not expected_labels:
        return Check("tracks_present", "SKIP", details="neither --tracks nor --track-config provided")
    if session_dict is None or not session_dict:
        return Check("tracks_present", "FAIL", details="sessionDictionary missing or empty")
    # Decode the first available session entry. Track names are identical
    # across per-region sessions (only the data: URL slices differ).
    sample_key = sorted(session_dict.keys())[0]
    session = decode_session_entry(session_dict[sample_key])
    if session is None:
        return Check("tracks_present", "FAIL", details="failed to decode/gunzip session entry")
    session_track_names = {t.get("name") for t in session.get("tracks", []) if t.get("name")}
    misses = [lab for lab in expected_labels if lab not in session_track_names]
    if misses:
        return Check(
            "tracks_present", "FAIL",
            observed=f"{len(expected_labels) - len(misses)}/{len(expected_labels)} found",
            expected=f"{len(expected_labels)}/{len(expected_labels)} found",
            details="missing: " + ", ".join(misses[:5]) + (" ..." if len(misses) > 5 else ""),
        )
    return Check(
        "tracks_present", "PASS",
        observed=f"{len(expected_labels)}/{len(expected_labels)} found",
    )


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def write_tsv(checks: list[Check], out: Path | None) -> None:
    lines = ["check\tstatus\tobserved\texpected\tdetails"]
    for c in checks:
        lines.append(f"{c.name}\t{c.status}\t{c.observed}\t{c.expected}\t{c.details}")
    text = "\n".join(lines) + "\n"
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    # Always also emit to stdout for piping / inspection.
    sys.stdout.write(text)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--html", required=True, help="path to the create_report HTML to verify")
    ap.add_argument("--sites", required=True, help="path to the sites BED that was passed to create_report")
    ap.add_argument("--tracks", nargs="*", default=[], help="track paths that were passed to create_report (--tracks mode)")
    ap.add_argument("--track-config", help="track config JSON that was passed to create_report (--track-config mode)")
    ap.add_argument("--min-size-mb", type=float, default=0.5, help="minimum acceptable HTML size in MB (default: 0.5)")
    ap.add_argument("--out", help="write the TSV report here in addition to stdout")
    ap.add_argument("--fail-on-fail", action="store_true", help="exit nonzero if any check is FAIL")
    args = ap.parse_args()

    html_path = Path(args.html)
    sites_path = Path(args.sites)
    out_path = Path(args.out) if args.out else None
    track_config = Path(args.track_config) if args.track_config else None

    checks: list[Check] = [check_html_exists(html_path)]

    # If the HTML doesn't exist, every downstream check would crash; mark them SKIP and bail.
    if checks[0].status == "FAIL":
        checks.append(Check("html_min_size", "SKIP", details="HTML missing"))
        checks.append(Check("region_count", "SKIP", details="HTML missing"))
        checks.append(Check("region_coords", "SKIP", details="HTML missing"))
        checks.append(Check("region_sessions", "SKIP", details="HTML missing"))
        checks.append(Check("tracks_present", "SKIP", details="HTML missing"))
        write_tsv(checks, out_path)
        if args.fail_on_fail:
            sys.exit(1)
        return

    checks.append(check_html_min_size(html_path, args.min_size_mb))

    html_text = html_path.read_text()
    table_json = parse_table_json(html_text)
    session_dict = parse_session_dictionary(html_text)
    bed_rows = load_sites_bed(sites_path)

    checks.append(check_region_count(bed_rows, table_json))
    checks.append(check_region_coords(bed_rows, table_json))
    checks.append(check_region_sessions(table_json, session_dict))
    checks.append(check_tracks_present(session_dict, expected_track_labels(args.tracks, track_config)))

    write_tsv(checks, out_path)

    if args.fail_on_fail and any(c.status == "FAIL" for c in checks):
        sys.exit(1)


if __name__ == "__main__":
    main()
