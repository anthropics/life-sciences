#!/usr/bin/env python3
"""Validate narrowPeak files from ENCODE ATAC-seq/DNase-seq accessibility aggregation.

Checks narrowPeak format, coordinate validity, typical accessibility peak sizes,
Tn5 bias artifacts, and reports summary statistics.

Usage:
    python validate_peaks.py input.narrowPeak [--blacklist hg38-blacklist.v2.bed] [--assay atac|dnase]
    python validate_peaks.py input.narrowPeak --assay atac --blacklist hg38-blacklist.v2.bed

Plain and gzipped (.gz) inputs and blacklists are both accepted.
"""

import argparse
import gzip
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

VALID_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
NARROWPEAK_COLS = 10

# Typical ATAC-seq peak widths: 100-500bp (nucleosome-free regions)
# Peaks >2kb are suspicious for accessibility data
ATAC_LARGE_THRESHOLD = 2000
DNASE_LARGE_THRESHOLD = 3000
ATAC_SMALL_THRESHOLD = 50  # Sub-nucleosomal Tn5 artifacts

MAX_COLUMN_ERRORS = 5
MAX_WARNINGS = 20


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate narrowPeak files from ATAC-seq/DNase-seq accessibility aggregation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python validate_peaks.py sample.narrowPeak\n"
            "  python validate_peaks.py sample.narrowPeak --assay atac\n"
            "  python validate_peaks.py sample.narrowPeak --blacklist hg38-blacklist.v2.bed\n"
        ),
    )
    parser.add_argument("input", type=Path, help="Input narrowPeak file")
    parser.add_argument(
        "--blacklist",
        type=Path,
        default=None,
        help="ENCODE blacklist BED file (e.g., hg38-blacklist.v2.bed)",
    )
    parser.add_argument(
        "--assay",
        choices=["atac", "dnase", "unknown"],
        default="unknown",
        help="Assay type for assay-specific checks. Default: unknown",
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


def validate_accessibility_peaks(input_path, blacklist_path, assay):
    errors = []
    warnings = []
    dropped_warnings = 0
    chrom_counts = Counter()
    peak_sizes = []
    signal_values = []
    start_positions = Counter()  # Track exact start positions for Tn5 pileup detection
    total_lines = 0
    # Every line excluded from the statistics below counts as malformed, so
    # total_lines == valid_peaks + bad_lines always holds.
    bad_lines = 0
    column_errors = 0
    blacklist_overlaps = 0
    chrm_peaks = 0
    large_peaks = 0
    tiny_peaks = 0

    large_threshold = ATAC_LARGE_THRESHOLD if assay == "atac" else DNASE_LARGE_THRESHOLD

    blacklist = None
    if blacklist_path:
        if not blacklist_path.exists():
            print(f"ERROR: Blacklist file not found: {blacklist_path}", file=sys.stderr)
            sys.exit(1)
        blacklist = load_blacklist(blacklist_path)

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

            # Accessibility peaks should always be narrowPeak
            if len(fields) != NARROWPEAK_COLS:
                bad_lines += 1
                column_errors += 1
                if column_errors <= MAX_COLUMN_ERRORS:
                    errors.append(
                        f"Line {line_num}: expected {NARROWPEAK_COLS} columns (narrowPeak), got {len(fields)}"
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
                if len(warnings) < MAX_WARNINGS:
                    warnings.append(
                        f"Line {line_num}: non-standard chromosome '{chrom}' (not in chr1-22, chrX, chrY, chrM)"
                    )
                else:
                    dropped_warnings += 1

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

            # signalValue, pValue and qValue are checked before anything is counted: a row with
            # an unusable value is malformed and stays out of every statistic
            scores = {}
            for col_idx, col_name in [(6, "signalValue"), (7, "pValue"), (8, "qValue")]:
                try:
                    scores[col_name] = float(fields[col_idx])
                except ValueError:
                    errors.append(f"Line {line_num}: invalid {col_name} in column {col_idx + 1}")
            if scores.get("signalValue", 0) < 0:
                errors.append(f"Line {line_num}: negative signalValue ({scores['signalValue']})")
            if len(scores) < 3 or scores["signalValue"] < 0:
                bad_lines += 1
                continue

            peak_size = end - start
            peak_sizes.append(peak_size)
            chrom_counts[chrom] += 1
            signal_values.append(scores["signalValue"])

            # Track start positions for Tn5 pileup detection (ATAC-specific)
            if assay == "atac":
                start_positions[(chrom, start)] += 1

            if chrom == "chrM":
                chrm_peaks += 1

            if peak_size > large_threshold:
                large_peaks += 1

            if peak_size < ATAC_SMALL_THRESHOLD:
                tiny_peaks += 1

            # Blacklist overlap check
            if blacklist and overlaps_blacklist(chrom, start, end, blacklist):
                blacklist_overlaps += 1

    if total_lines == 0:
        print(f"ERROR: no data rows in {input_path} (only comments, headers or blank lines)", file=sys.stderr)
        sys.exit(1)

    valid_peaks = total_lines - bad_lines

    # --- Tn5 Pileup Detection (ATAC-specific) ---
    tn5_pileup_count = 0
    if assay == "atac" and start_positions:
        tn5_pileup_count = sum(1 for count in start_positions.values() if count >= 5)

    # --- Report Statistics ---
    assay_label = assay.upper() if assay != "unknown" else "Accessibility"
    print(f"=== {assay_label} NarrowPeak Validation Report ===")
    print(f"File: {input_path}")
    print(f"Assay: {assay}")
    print()

    print("--- Summary ---")
    print(f"Data lines: {total_lines:,}")
    print(f"Valid peaks: {valid_peaks:,}")
    print(f"Malformed lines: {bad_lines}")
    if blacklist_path:
        print(f"Blacklist overlaps: {blacklist_overlaps:,} ({100 * blacklist_overlaps / max(valid_peaks, 1):.1f}%)")
    print()

    if peak_sizes:
        sorted_sizes = sorted(peak_sizes)
        n = len(sorted_sizes)
        size_q1, size_median, size_q3 = quartiles(sorted_sizes)

        # Count peaks in expected accessibility range (100-500bp)
        in_range = sum(1 for s in peak_sizes if 100 <= s <= 500)
        in_range_pct = 100 * in_range / max(n, 1)

        print("--- Peak Size Distribution ---")
        print(f"Min:    {sorted_sizes[0]:,} bp")
        print(f"25th:   {size_q1:,.1f} bp")
        print(f"Median: {size_median:,.1f} bp")
        print(f"75th:   {size_q3:,.1f} bp")
        print(f"Max:    {sorted_sizes[-1]:,} bp")
        print(f"In typical range (100-500bp): {in_range:,} ({in_range_pct:.1f}%)")
        print()

    if signal_values:
        sorted_sig = sorted(signal_values)
        sig_q1, sig_median, sig_q3 = quartiles(sorted_sig)
        print("--- SignalValue Distribution ---")
        print(f"Min:    {sorted_sig[0]:.2f}")
        print(f"25th:   {sig_q1:.2f}")
        print(f"Median: {sig_median:.2f}")
        print(f"75th:   {sig_q3:.2f}")
        print(f"Max:    {sorted_sig[-1]:.2f}")
        print()

    print("--- Chromosome Distribution ---")
    for chrom in sorted(chrom_counts.keys(), key=lambda c: (len(c), c)):
        count = chrom_counts[chrom]
        pct = 100 * count / max(valid_peaks, 1)
        print(f"  {chrom:<6} {count:>8,}  ({pct:5.1f}%)")
    print()

    # --- Warnings ---
    if chrm_peaks > 0:
        msg = (
            f"WARNING: {chrm_peaks:,} peaks on chrM. "
            f"Mitochondrial peaks are common ATAC-seq artifacts (high mito read fraction)."
        )
        print(msg, file=sys.stderr)

    if large_peaks > 0:
        msg = (
            f"WARNING: {large_peaks:,} peaks exceed {large_threshold}bp. "
            f"Accessibility peaks are typically 100-500bp. "
            f"Large peaks may indicate artifacts or broad-mark contamination."
        )
        print(msg, file=sys.stderr)

    if tiny_peaks > 0 and assay == "atac":
        msg = (
            f"WARNING: {tiny_peaks:,} peaks are <{ATAC_SMALL_THRESHOLD}bp. "
            f"Very narrow peaks in ATAC-seq can be Tn5 insertion artifacts."
        )
        print(msg, file=sys.stderr)

    if tn5_pileup_count > 0:
        msg = (
            f"WARNING: {tn5_pileup_count:,} genomic positions have 5+ peaks "
            f"sharing the exact same start coordinate. This may indicate Tn5 "
            f"insertion bias (positional pileup artifact)."
        )
        print(msg, file=sys.stderr)

    if blacklist_overlaps > 0:
        msg = (
            f"WARNING: {blacklist_overlaps:,} peaks overlap ENCODE blacklist regions. Remove these before aggregation."
        )
        print(msg, file=sys.stderr)

    for w in warnings:
        print(w, file=sys.stderr)
    if dropped_warnings:
        print(f"  ... and {dropped_warnings} more warning(s) suppressed", file=sys.stderr)

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
        print("\nRESULT: PASS -- file is valid narrowPeak")

    return 1 if has_errors else 0


if __name__ == "__main__":
    args = parse_args()
    exit_code = validate_accessibility_peaks(args.input, args.blacklist, args.assay)
    sys.exit(exit_code)
