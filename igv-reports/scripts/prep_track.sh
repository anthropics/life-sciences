#!/usr/bin/env bash
# prep_track.sh — convert a plain-gzip GFF3/GTF/BED.gz into a properly
# bgzipped + tabix-indexed track that igv-reports can load.
#
# Author: Samuel Ahuno
# Why: igv-reports parses tracks by extension and needs bgzip+tabix.
# Plain gzip with `.gz` extension trips it with a UnicodeDecodeError or
# silently fails. Tabix indexing additionally requires position-sorted
# records within each chromosome, which gencode/many-other distributions
# do not guarantee — they interleave records by feature type.
#
# Pipeline: backup -> gunzip -> sort by chr+pos (preserving header) ->
# bgzip in place -> tabix -p <gff|gtf|bed>.
#
# Usage:
#   prep_track.sh <track.gff3.gz | track.gtf.gz | track.bed.gz>
#
# Output:
#   <input>                              (replaced with new bgzip)
#   <input>.tbi                          (new tabix index)
#   <input>.bak.original_gzip            (backup of the original .gz)

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <track.gff3.gz | track.gtf.gz | track.bed.gz>" >&2
    exit 2
fi

INPUT=$1
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: file not found: $INPUT" >&2
    exit 2
fi

# Detect format by suffix.
case "$INPUT" in
    *.gff3.gz|*.gff.gz)  FMT=gff ;;
    *.gtf.gz)            FMT=gff ;;   # tabix preset for GTF is named "gff"
    *.bed.gz|*.bedgraph.gz) FMT=bed ;;
    *) echo "ERROR: unsupported extension: $INPUT (need .gff3.gz, .gtf.gz, .bed.gz, .bedgraph.gz)" >&2; exit 2 ;;
esac

# Need bgzip / tabix / sort / gunzip.
for tool in bgzip tabix sort gunzip awk file; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool not on PATH. Activate the snakemake conda env first." >&2
        exit 2
    fi
done

# Detect if already bgzip — skip the whole conversion if it is and just
# rebuild the index.
if file "$INPUT" | grep -q "extra field"; then
    echo "[$(date '+%F %T')] $INPUT is already bgzip; rebuilding tabix index only."
    rm -f "${INPUT}.tbi"
    tabix -p "$FMT" "$INPUT"
    echo "[$(date '+%F %T')] DONE: ${INPUT}.tbi"
    exit 0
fi

BACKUP="${INPUT}.bak.original_gzip"
if [[ -f "$BACKUP" ]]; then
    echo "[$(date '+%F %T')] backup already exists: $BACKUP — refusing to overwrite. Move it aside and rerun if you want a fresh backup."
else
    cp -p "$INPUT" "$BACKUP"
    echo "[$(date '+%F %T')] backed up to $BACKUP"
fi

# Decompress to a sibling temp.
TMP="${INPUT%.gz}.unsorted.tmp"
SORTED="${INPUT%.gz}.sorted.tmp"
gunzip -c "$INPUT" > "$TMP"
echo "[$(date '+%F %T')] decompressed to $TMP ($(stat -c %s "$TMP") bytes)"

# Sort: preserve any leading # header lines, sort body by chr (column 1)
# then numeric pos (column 4 for GFF/GTF; column 2 for BED).
case "$FMT" in
    gff)  POS_COL=4 ;;
    bed)  POS_COL=2 ;;
esac

(grep '^#' "$TMP" || true) > "$SORTED"
grep -v '^#' "$TMP" \
    | sort -k1,1 -k${POS_COL},${POS_COL}n -S 2G --parallel=4 \
    >> "$SORTED"
echo "[$(date '+%F %T')] sorted by chr,pos (col $POS_COL) into $SORTED"

# bgzip into place (replaces the original) and index.
mv "$SORTED" "${INPUT%.gz}"
rm -f "$TMP"
bgzip -@ 4 "${INPUT%.gz}"
echo "[$(date '+%F %T')] bgzipped: $INPUT ($(stat -c %s "$INPUT") bytes)"

rm -f "${INPUT}.tbi"
tabix -p "$FMT" "$INPUT"
echo "[$(date '+%F %T')] indexed: ${INPUT}.tbi ($(stat -c %s "${INPUT}.tbi") bytes)"

# Sanity check: pull the first contig's first 100 kb and confirm tabix returns rows.
FIRST_CONTIG=$(zcat "$INPUT" | awk '$1!~/^#/ {print $1; exit}')
if [[ -n "$FIRST_CONTIG" ]]; then
    N=$(tabix "$INPUT" "${FIRST_CONTIG}:1-100000" | wc -l)
    echo "[$(date '+%F %T')] sanity: ${FIRST_CONTIG}:1-100000 returns $N row(s)"
fi

echo "[$(date '+%F %T')] DONE — track ready for igv-reports. Original preserved at $BACKUP"
