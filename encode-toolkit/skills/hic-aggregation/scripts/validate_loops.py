#!/usr/bin/env python3
"""Validate BEDPE loop files from ENCODE Hi-C aggregation.

Checks BEDPE format, anchor validity, resolution consistency, canonical ordering,
cis/trans classification, and reports summary statistics.

Usage:
    python validate_loops.py input.bedpe [--min-distance 20000] [--expected-resolution 10000]
    python validate_loops.py input.bedpe --expected-resolution 5000

Plain and gzipped (.gz) inputs are both accepted.
"""

import argparse
import gzip
import statistics
import sys
from collections import Counter
from pathlib import Path

VALID_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
MIN_BEDPE_COLS = 6

MAX_COLUMN_ERRORS = 5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate BEDPE loop files from ENCODE Hi-C aggregation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python validate_loops.py loops.bedpe\n"
            "  python validate_loops.py loops.bedpe --expected-resolution 10000\n"
            "  python validate_loops.py loops.bedpe --min-distance 20000\n"
        ),
    )
    parser.add_argument("input", type=Path, help="Input BEDPE loop file")
    parser.add_argument(
        "--min-distance",
        type=int,
        default=20000,
        help="Minimum anchor-to-anchor distance for cis loops (bp). Default: 20000",
    )
    parser.add_argument(
        "--expected-resolution",
        type=int,
        default=None,
        help="Expected resolution in bp (e.g., 5000, 10000, 25000). Checks anchor size consistency.",
    )
    return parser.parse_args()


def open_text(path):
    """Open a plain or gzipped text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def fmt_distance(bp):
    """Render a genomic distance, in bp below 1kb so sub-kb values do not read as 0kb."""
    return f"{bp}bp" if bp < 1000 else f"{bp // 1000}kb"


def quartiles(values):
    """25th percentile, median and 75th percentile (interpolated)."""
    median = statistics.median(values)
    if len(values) < 2:
        return values[0], median, values[0]
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive")
    return q1, median, q3


def validate_loops(input_path, min_distance, expected_resolution):
    errors = []
    warnings = []
    chrom_counts = Counter()
    anchor1_sizes = []
    anchor2_sizes = []
    loop_distances = []
    total_lines = 0
    # Every line excluded from the statistics below counts as malformed, so
    # total_lines == valid_loops + bad_lines always holds.
    bad_lines = 0
    column_errors = 0
    cis_loops = 0
    trans_loops = 0
    short_range = 0
    non_canonical = 0

    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    with open_text(input_path) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            line = line.strip()
            if not line:
                continue

            total_lines += 1
            fields = line.split("\t")

            # BEDPE requires at least 6 columns
            if len(fields) < MIN_BEDPE_COLS:
                bad_lines += 1
                column_errors += 1
                if column_errors <= MAX_COLUMN_ERRORS:
                    errors.append(f"Line {line_num}: expected >= {MIN_BEDPE_COLS} columns (BEDPE), got {len(fields)}")
                elif column_errors == MAX_COLUMN_ERRORS + 1:
                    errors.append("... suppressing further column-count errors")
                continue

            chr1 = fields[0]
            chr2 = fields[3]

            # Chromosome validation for both anchors: a row with either anchor on an
            # unrecognized chromosome is one malformed line, and is skipped entirely.
            invalid_anchors = [
                (label, value)
                for label, value in (("anchor1", chr1), ("anchor2", chr2))
                if value not in VALID_CHROMS and not value.startswith("chr")
            ]
            if invalid_anchors:
                for label, value in invalid_anchors:
                    errors.append(f"Line {line_num}: invalid {label} chromosome '{value}'")
                bad_lines += 1
                continue

            # Coordinate validation
            try:
                start1 = int(fields[1])
                end1 = int(fields[2])
                start2 = int(fields[4])
                end2 = int(fields[5])
            except ValueError:
                errors.append(f"Line {line_num}: non-integer coordinates in anchor fields")
                bad_lines += 1
                continue

            # Non-negative coordinates, and start < end for each anchor. A row that breaks either
            # rule is malformed: count it once and keep it out of the statistics.
            impossible_anchor = False
            for label, val in [("start1", start1), ("end1", end1), ("start2", start2), ("end2", end2)]:
                if val < 0:
                    errors.append(f"Line {line_num}: negative {label} ({val})")
                    impossible_anchor = True
            if start1 >= end1:
                errors.append(f"Line {line_num}: anchor1 start ({start1}) >= end ({end1})")
                impossible_anchor = True
            if start2 >= end2:
                errors.append(f"Line {line_num}: anchor2 start ({start2}) >= end ({end2})")
                impossible_anchor = True
            if impossible_anchor:
                bad_lines += 1
                continue

            a1_size = end1 - start1
            a2_size = end2 - start2
            anchor1_sizes.append(a1_size)
            anchor2_sizes.append(a2_size)

            # Cis vs trans classification
            if chr1 == chr2:
                cis_loops += 1
                chrom_counts[chr1] += 1

                # Distance check (midpoint-to-midpoint)
                mid1 = (start1 + end1) // 2
                mid2 = (start2 + end2) // 2
                distance = abs(mid2 - mid1)
                loop_distances.append(distance)

                if distance < min_distance:
                    short_range += 1

                # Canonical ordering: anchor1.start < anchor2.start
                if start1 > start2:
                    non_canonical += 1
            else:
                trans_loops += 1

            # Expected resolution check
            if expected_resolution is not None:
                if a1_size != expected_resolution:
                    if len(warnings) < 10:
                        warnings.append(f"Line {line_num}: anchor1 size {a1_size} != expected {expected_resolution}")
                if a2_size != expected_resolution:
                    if len(warnings) < 10:
                        warnings.append(f"Line {line_num}: anchor2 size {a2_size} != expected {expected_resolution}")

    if total_lines == 0:
        print(f"ERROR: no data rows in {input_path} (only comments, headers or blank lines)", file=sys.stderr)
        sys.exit(1)

    valid_loops = total_lines - bad_lines

    # --- Detect resolution from anchor sizes ---
    detected_resolution = None
    if anchor1_sizes:
        all_sizes = anchor1_sizes + anchor2_sizes
        size_counts = Counter(all_sizes)
        most_common_size, most_common_count = size_counts.most_common(1)[0]
        size_pct = 100 * most_common_count / len(all_sizes)
        if size_pct > 50:
            detected_resolution = most_common_size

    # --- Report Statistics ---
    print("=== Hi-C BEDPE Loop Validation Report ===")
    print(f"File: {input_path}")
    print()

    print("--- Summary ---")
    print(f"Data lines: {total_lines:,}")
    print(f"Valid loops: {valid_loops:,}")
    print(f"Malformed lines: {bad_lines}")
    print(f"Cis loops (same chromosome): {cis_loops:,} ({100 * cis_loops / max(valid_loops, 1):.1f}%)")
    print(f"Trans loops (inter-chromosomal): {trans_loops:,} ({100 * trans_loops / max(valid_loops, 1):.1f}%)")
    print(f"Non-canonical ordering (anchor1 > anchor2): {non_canonical:,}")
    print(f"Short-range loops (<{fmt_distance(min_distance)}): {short_range:,}")
    print()

    if detected_resolution:
        print("--- Resolution ---")
        print(f"Detected resolution: {detected_resolution:,} bp ({fmt_distance(detected_resolution)})")
        print(f"Anchor size consistency: {size_pct:.1f}% of anchors match detected resolution")
    elif anchor1_sizes:
        print("--- Resolution ---")
        print("WARNING: No consistent resolution detected. Anchor sizes vary.")
        top3 = size_counts.most_common(3)
        for sz, cnt in top3:
            print(f"  {sz:>8,} bp: {cnt:,} anchors ({100 * cnt / len(all_sizes):.1f}%)")
    print()

    if loop_distances:
        sorted_dist = sorted(loop_distances)
        n = len(sorted_dist)
        dist_q1, dist_median, dist_q3 = quartiles(sorted_dist)
        print("--- Loop Distance Distribution (cis only) ---")
        print(f"Min:    {sorted_dist[0]:>12,} bp ({fmt_distance(sorted_dist[0])})")
        print(f"25th:   {dist_q1:>12,.1f} bp ({fmt_distance(int(dist_q1))})")
        print(f"Median: {dist_median:>12,.1f} bp ({fmt_distance(int(dist_median))})")
        print(f"75th:   {dist_q3:>12,.1f} bp ({fmt_distance(int(dist_q3))})")
        print(f"Max:    {sorted_dist[-1]:>12,} bp ({fmt_distance(sorted_dist[-1])})")

        # Distance buckets
        buckets = [
            ("< 100kb", 0, 100_000),
            ("100kb - 500kb", 100_000, 500_000),
            ("500kb - 1Mb", 500_000, 1_000_000),
            ("1Mb - 5Mb", 1_000_000, 5_000_000),
            (">= 5Mb", 5_000_000, float("inf")),
        ]
        print("\n  Distance buckets:")
        for label, lo, hi in buckets:
            count = sum(1 for d in loop_distances if lo <= d < hi)
            print(f"    {label:<15} {count:>8,}  ({100 * count / n:.1f}%)")
        print()

    if cis_loops > 0:
        print("--- Chromosome Distribution (cis loops) ---")
        for chrom in sorted(chrom_counts.keys(), key=lambda c: (len(c), c)):
            count = chrom_counts[chrom]
            pct = 100 * count / max(cis_loops, 1)
            print(f"  {chrom:<6} {count:>8,}  ({pct:5.1f}%)")
        print()

    # --- Warnings ---
    if trans_loops > valid_loops * 0.05:
        msg = (
            f"WARNING: {100 * trans_loops / max(valid_loops, 1):.1f}% of loops are trans "
            f"(inter-chromosomal). Expect <5% for typical Hi-C data."
        )
        print(msg, file=sys.stderr)

    if short_range > 0:
        msg = (
            f"WARNING: {short_range:,} loops have anchors <{fmt_distance(min_distance)} apart. "
            f"These may be self-ligation artifacts."
        )
        print(msg, file=sys.stderr)

    if non_canonical > 0:
        msg = (
            f"WARNING: {non_canonical:,} loops have non-canonical ordering "
            f"(anchor1.start > anchor2.start). Canonicalize before merging."
        )
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
        print("\nRESULT: PASS -- file is valid BEDPE")

    return 1 if has_errors else 0


if __name__ == "__main__":
    args = parse_args()
    exit_code = validate_loops(args.input, args.min_distance, args.expected_resolution)
    sys.exit(exit_code)
