#!/usr/bin/env python3
"""Validate bedMethyl files from ENCODE WGBS methylation aggregation.

Checks bedMethyl format, methylation value ranges, strand validity, coverage
values, and reports summary statistics with warnings for low-coverage sites.

Usage:
    python validate_methylation.py input.bedMethyl [--min-coverage 5] [--blacklist hg38-blacklist.v2.bed]
    python validate_methylation.py input.bed --min-coverage 10
    python validate_methylation.py input.bedMethyl --scale fraction

Plain and gzipped (.gz) inputs and blacklists are both accepted.
"""

import argparse
import gzip
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

VALID_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
VALID_STRANDS = {"+", "-", "."}

# ENCODE bedMethyl format has 11 columns:
# chr start end name score strand thickStart thickEnd color coverage percentMethylated
BEDMETHYL_COLS = 11

# Alternative minimal format: chr start end name score strand coverage methylation%
MINIMAL_COLS = 8

# layout -> (coverage column, methylation column, how the column count is described), 0-indexed
LAYOUTS = {
    "encode_bedmethyl": (9, 10, "11+"),
    "minimal": (6, 7, "8"),
}

MAX_COLUMN_ERRORS = 5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate bedMethyl files from ENCODE WGBS methylation aggregation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python validate_methylation.py sample.bedMethyl\n"
            "  python validate_methylation.py sample.bedMethyl --min-coverage 10\n"
            "  python validate_methylation.py sample.bed --blacklist hg38-blacklist.v2.bed\n"
        ),
    )
    parser.add_argument("input", type=Path, help="Input bedMethyl file")
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=5,
        help="Minimum coverage threshold to flag low-coverage CpGs. Default: 5",
    )
    parser.add_argument(
        "--blacklist",
        type=Path,
        default=None,
        help="ENCODE blacklist BED file (e.g., hg38-blacklist.v2.bed)",
    )
    parser.add_argument(
        "--scale",
        choices=["auto", "percent", "fraction"],
        default="auto",
        help=(
            "Scale of the methylation column, decided once for the whole file. "
            "auto: percent for the ENCODE 11-column layout, and for the minimal "
            "layout fraction only when the file's maximum value is <= 1. Default: auto"
        ),
    )
    return parser.parse_args()


def open_text(path):
    """Open a plain or gzipped text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def quartiles(values):
    """25th percentile, median and 75th percentile (interpolated)."""
    median = statistics.median(values)
    if len(values) < 2:
        return values[0], median, values[0]
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return q1, median, q3


def load_blacklist(path):
    """Load blacklist regions as a dict of chrom -> list of (start, end)."""
    regions = defaultdict(list)
    with open_text(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                regions[parts[0]].append((int(parts[1]), int(parts[2])))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def overlaps_blacklist(chrom, start, end, blacklist):
    """Check if a region overlaps any blacklist interval."""
    if chrom not in blacklist:
        return False
    for bl_start, bl_end in blacklist[chrom]:
        if bl_start >= end:
            break
        if bl_end > start:
            return True
    return False


def detect_format(n_cols):
    """ENCODE bedMethyl (11 or more columns), minimal (exactly 8), or None for anything else."""
    if n_cols >= BEDMETHYL_COLS:
        return "encode_bedmethyl"
    if n_cols == MINIMAL_COLS:
        return "minimal"
    return None


def layout_matches(detected_format, n_cols):
    """Check a later line against the layout detected from the first data line."""
    if detected_format == "encode_bedmethyl":
        return n_cols >= BEDMETHYL_COLS
    return n_cols == MINIMAL_COLS


def decide_scale(requested, detected_format, values):
    """Decide once per file whether the methylation column holds percentages or fractions."""
    if requested != "auto":
        return requested
    # ENCODE bedMethyl column 11 is a percentage (0-100) by specification.
    if detected_format == "encode_bedmethyl":
        return "percent"
    return "fraction" if values and max(values) <= 1.0 else "percent"


def validate_methylation(input_path, min_coverage, blacklist_path, requested_scale):
    errors = []
    warnings = []
    chrom_counts = Counter()
    strand_counts = Counter()
    coverage_values = []
    records = []  # rows whose fields all parsed: (line, chrom, strand, coverage, methylation, in blacklist)
    total_lines = 0
    # Every line excluded from the statistics below counts as malformed, so
    # total_lines == valid_records + bad_lines always holds.
    bad_lines = 0
    column_errors = 0
    low_coverage = 0
    blacklist_overlaps = 0

    blacklist = None
    if blacklist_path:
        if not blacklist_path.exists():
            print(f"ERROR: Blacklist file not found: {blacklist_path}", file=sys.stderr)
            sys.exit(1)
        blacklist = load_blacklist(blacklist_path)

    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    detected_format = None
    expected_desc = None
    cov_col = None  # 0-indexed column for coverage
    meth_col = None  # 0-indexed column for methylation

    with open_text(input_path) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            line = line.strip()
            if not line:
                continue

            total_lines += 1
            fields = line.split("\t")

            # The first data line fixes the layout for the whole file
            if detected_format is None:
                detected_format = detect_format(len(fields))
                if detected_format is None:
                    print(
                        f"ERROR: Could not detect bedMethyl format. "
                        f"Expected 8 columns (minimal) or >= 11 columns (ENCODE bedMethyl), got {len(fields)}",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                cov_col, meth_col, expected_desc = LAYOUTS[detected_format]

            if not layout_matches(detected_format, len(fields)):
                bad_lines += 1
                column_errors += 1
                if column_errors <= MAX_COLUMN_ERRORS:
                    errors.append(
                        f"Line {line_num}: expected {expected_desc} columns ({detected_format}), got {len(fields)}"
                    )
                elif column_errors == MAX_COLUMN_ERRORS + 1:
                    errors.append("... suppressing further column-count errors")
                continue

            chrom = fields[0]
            if chrom not in VALID_CHROMS:
                if not chrom.startswith("chr"):
                    errors.append(f"Line {line_num}: invalid chromosome '{chrom}'")
                    bad_lines += 1
                    continue

            # Coordinate validation
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                errors.append(f"Line {line_num}: non-integer coordinates")
                bad_lines += 1
                continue

            if start < 0:
                errors.append(f"Line {line_num}: negative start coordinate ({start})")
            if end < 0:
                errors.append(f"Line {line_num}: negative end coordinate ({end})")
            if start >= end:
                errors.append(f"Line {line_num}: start ({start}) >= end ({end})")
            # an impossible interval is malformed: count it once and keep it out of the statistics
            if start < 0 or end < 0 or start >= end:
                bad_lines += 1
                continue

            # Every field is checked before anything is counted: a row with an unusable strand,
            # coverage or methylation value is malformed and stays out of every statistic.
            row_ok = True
            strand = fields[5]
            if strand not in VALID_STRANDS:
                errors.append(f"Line {line_num}: invalid strand '{strand}' (expected +, -, or .)")
                row_ok = False

            try:
                coverage = int(float(fields[cov_col]))  # some tools write coverage as a float
            except ValueError:
                errors.append(f"Line {line_num}: invalid coverage in column {cov_col + 1}")
                row_ok = False
            else:
                if coverage < 0:
                    errors.append(f"Line {line_num}: negative coverage ({coverage})")
                    row_ok = False

            try:
                meth = float(fields[meth_col])
            except ValueError:
                errors.append(f"Line {line_num}: invalid methylation value in column {meth_col + 1}")
                row_ok = False
            else:
                if meth < 0:
                    errors.append(f"Line {line_num}: negative methylation value ({meth})")
                    row_ok = False

            if not row_ok:
                bad_lines += 1
                continue

            # The methylation scale is known only after the whole file is read, so the row is
            # kept here and counted below, once its value is known to be in range.
            in_blacklist = bool(blacklist) and overlaps_blacklist(chrom, start, end, blacklist)
            records.append((line_num, chrom, strand, coverage, meth, in_blacklist))

    if total_lines == 0:
        print(f"ERROR: no data rows in {input_path} (only comments, headers or blank lines)", file=sys.stderr)
        sys.exit(1)

    # --- Apply one scale to every methylation value, then count the rows that are in range ---
    scale = decide_scale(requested_scale, detected_format, [record[4] for record in records])
    upper, factor = (1.0, 100) if scale == "fraction" else (100, 1)
    methylation_values = []
    for line_num, chrom, strand, coverage, meth, in_blacklist in records:
        if meth > upper:
            if scale == "fraction":
                errors.append(f"Line {line_num}: methylation value {meth} > 1 with --scale fraction (expected 0-1).")
            else:
                errors.append(f"Line {line_num}: methylation value > 100 ({meth}). Expected 0-100 (percentage).")
            bad_lines += 1
            continue
        chrom_counts[chrom] += 1
        strand_counts[strand] += 1
        coverage_values.append(coverage)
        if coverage < min_coverage:
            low_coverage += 1
        if in_blacklist:
            blacklist_overlaps += 1
        methylation_values.append(meth * factor)

    valid_records = total_lines - bad_lines

    # --- Report Statistics ---
    print("=== bedMethyl Validation Report ===")
    print(f"File: {input_path}")
    print(f"Detected format: {detected_format} ({expected_desc} columns)")
    print()

    print("--- Summary ---")
    print(f"Data lines: {total_lines:,}")
    print(f"Valid CpGs: {valid_records:,}")
    print(f"Malformed lines: {bad_lines}")
    print(f"Low-coverage CpGs (<{min_coverage}x): {low_coverage:,} ({100 * low_coverage / max(valid_records, 1):.1f}%)")
    if blacklist_path:
        print(f"Blacklist overlaps: {blacklist_overlaps:,} ({100 * blacklist_overlaps / max(valid_records, 1):.1f}%)")
    print()

    scale_source = "auto-detected" if requested_scale == "auto" else f"--scale {requested_scale}"
    scale_label = "fraction (0-1), reported as percent" if scale == "fraction" else "percent (0-100)"
    print(f"Methylation scale: {scale_label} [{scale_source}]")
    print()

    if coverage_values:
        sorted_cov = sorted(coverage_values)
        n = len(sorted_cov)
        cov_q1, cov_median, cov_q3 = quartiles(sorted_cov)
        print("--- Coverage Distribution ---")
        print(f"Min:    {sorted_cov[0]:>6}")
        print(f"25th:   {cov_q1:>6.1f}")
        print(f"Median: {cov_median:>6.1f}")
        print(f"75th:   {cov_q3:>6.1f}")
        print(f"Max:    {sorted_cov[-1]:>6}")

        # Coverage buckets
        buckets = [
            ("<3x", 0, 3),
            ("3-5x", 3, 5),
            ("5-10x", 5, 10),
            ("10-20x", 10, 20),
            ("20-50x", 20, 50),
            (">=50x", 50, float("inf")),
        ]
        print("\n  Coverage buckets:")
        for label, lo, hi in buckets:
            count = sum(1 for c in coverage_values if lo <= c < hi)
            print(f"    {label:<8} {count:>10,}  ({100 * count / n:.1f}%)")
        print()

    if methylation_values:
        sorted_meth = sorted(methylation_values)
        n_m = len(sorted_meth)
        meth_q1, meth_median, meth_q3 = quartiles(sorted_meth)
        print("--- Methylation Distribution (as %) ---")
        print(f"Min:    {sorted_meth[0]:>6.1f}%")
        print(f"25th:   {meth_q1:>6.1f}%")
        print(f"Median: {meth_median:>6.1f}%")
        print(f"75th:   {meth_q3:>6.1f}%")
        print(f"Max:    {sorted_meth[-1]:>6.1f}%")

        # Methylation state buckets
        buckets = [
            ("Unmethylated (0-10%)", 0, 10),
            ("Low (10-30%)", 10, 30),
            ("Intermediate (30-70%)", 30, 70),
            ("High (70-90%)", 70, 90),
            ("Methylated (90-100%)", 90, 100.01),
        ]
        print("\n  Methylation state distribution:")
        for label, lo, hi in buckets:
            count = sum(1 for m in methylation_values if lo <= m < hi)
            print(f"    {label:<30} {count:>10,}  ({100 * count / n_m:.1f}%)")
        print()

    print("--- Strand Breakdown ---")
    for strand in ["+", "-", "."]:
        count = strand_counts.get(strand, 0)
        pct = 100 * count / max(valid_records, 1)
        print(f"  {strand:<3} {count:>10,}  ({pct:5.1f}%)")
    print()

    print("--- Chromosome Distribution ---")
    for chrom in sorted(chrom_counts.keys(), key=lambda c: (len(c), c)):
        count = chrom_counts[chrom]
        pct = 100 * count / max(valid_records, 1)
        print(f"  {chrom:<6} {count:>10,}  ({pct:5.1f}%)")
    print()

    # --- Warnings ---
    if low_coverage > valid_records * 0.3:
        msg = (
            f"WARNING: {100 * low_coverage / max(valid_records, 1):.0f}% of CpGs have "
            f"coverage <{min_coverage}x. Consider filtering these for reliable "
            f"methylation estimates."
        )
        print(msg, file=sys.stderr)

    if valid_records and strand_counts.get(".", 0) == valid_records:
        msg = (
            "INFO: All CpGs have strand '.'. This file may already be strand-merged. "
            "Skip the strand-merge step in aggregation."
        )
        print(msg, file=sys.stderr)
    elif strand_counts.get("+", 0) > 0 and strand_counts.get("-", 0) > 0:
        plus_count = strand_counts.get("+", 0)
        minus_count = strand_counts.get("-", 0)
        ratio = plus_count / max(minus_count, 1)
        if 0.8 <= ratio <= 1.2:
            msg = (
                f"INFO: Both strands present ({plus_count:,} forward, {minus_count:,} reverse). "
                f"Consider strand-merging for increased per-CpG coverage."
            )
            print(msg, file=sys.stderr)

    if blacklist_overlaps > 0:
        msg = f"WARNING: {blacklist_overlaps:,} CpGs overlap ENCODE blacklist regions. Remove these before aggregation."
        print(msg, file=sys.stderr)

    for w in warnings[:20]:
        print(w, file=sys.stderr)

    # --- Errors ---
    if errors:
        print(f"\n--- Errors ({len(errors)}) ---", file=sys.stderr)
        for e in errors[:50]:
            print(f"  {e}", file=sys.stderr)
        if len(errors) > 50:
            print(f"  ... and {len(errors) - 50} more errors", file=sys.stderr)

    has_errors = len(errors) > 0
    if has_errors:
        print(f"\nRESULT: FAIL -- {len(errors)} error(s) found", file=sys.stderr)
    else:
        print(f"\nRESULT: PASS -- file is valid {detected_format}")

    return 1 if has_errors else 0


if __name__ == "__main__":
    args = parse_args()
    exit_code = validate_methylation(args.input, args.min_coverage, args.blacklist, args.scale)
    sys.exit(exit_code)
