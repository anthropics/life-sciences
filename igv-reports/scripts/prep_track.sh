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
#   prep_track.sh <input.gz> --out <sibling.bgz.gz>
#
# In-place mode (default):
#   <input>                              (replaced with new bgzip)
#   <input>.tbi                          (new tabix index)
#   <input>.bak.original_gzip            (backup of the original .gz)
#
# Sibling mode (--out PATH; non-destructive):
#   <input>                              (unchanged)
#   <out>                                (new bgzip — same extension family as input)
#   <out>.tbi                            (new tabix index)
#   (no backup created — original is left as-is)

set -euo pipefail

INPUT=""
OUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out)
            [[ $# -lt 2 ]] && { echo "ERROR: --out requires a path" >&2; exit 2; }
            OUT=$2; shift 2 ;;
        --out=*)
            OUT=${1#--out=}; shift ;;
        -h|--help)
            sed -n '2,28p' "$0" >&2; exit 0 ;;
        --)
            shift; break ;;
        -*)
            echo "ERROR: unknown flag: $1" >&2; exit 2 ;;
        *)
            if [[ -z "$INPUT" ]]; then INPUT=$1
            else echo "ERROR: unexpected positional arg: $1" >&2; exit 2
            fi
            shift ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "Usage: $0 <track.gff3.gz | track.gtf.gz | track.bed.gz> [--out <sibling.gz>]" >&2
    exit 2
fi
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: file not found: $INPUT" >&2
    exit 2
fi
if [[ -n "$OUT" && -e "$OUT" ]]; then
    echo "ERROR: --out target already exists: $OUT — refusing to overwrite. Move it aside and rerun." >&2
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
        echo "ERROR: $tool not on PATH. Install htslib and add bgzip/tabix to PATH first." >&2
        exit 2
    fi
done

# Resolve where the final bgzip + .tbi will land. In sibling mode we never
# touch the original. In in-place mode the target IS the original, with a
# backup taken first.
if [[ -n "$OUT" ]]; then
    TARGET=$OUT
    mkdir -p "$(dirname "$TARGET")"
else
    TARGET=$INPUT
fi

# Detect if already bgzip — skip the whole conversion if it is and just
# rebuild the index. (In sibling mode this means: copy + index, leaving
# the original untouched.)
if file "$INPUT" | grep -q "extra field"; then
    if [[ "$TARGET" != "$INPUT" ]]; then
        cp -p "$INPUT" "$TARGET"
        echo "[$(date '+%F %T')] $INPUT already bgzip; copied to $TARGET, rebuilding tabix index."
    else
        echo "[$(date '+%F %T')] $INPUT is already bgzip; rebuilding tabix index only."
    fi
    rm -f "${TARGET}.tbi"
    tabix -p "$FMT" "$TARGET"
    echo "[$(date '+%F %T')] DONE: ${TARGET}.tbi"
    exit 0
fi

# In-place mode: take a backup of the original. In sibling mode no backup is
# needed since the original is never modified.
if [[ "$TARGET" == "$INPUT" ]]; then
    BACKUP="${INPUT}.bak.original_gzip"
    if [[ -f "$BACKUP" ]]; then
        echo "[$(date '+%F %T')] backup already exists: $BACKUP — refusing to overwrite. Move it aside and rerun if you want a fresh backup."
    else
        cp -p "$INPUT" "$BACKUP"
        echo "[$(date '+%F %T')] backed up to $BACKUP"
    fi
fi

# Decompress to a sibling-of-INPUT temp (always, regardless of target).
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

# bgzip and index. Sibling mode: SORTED -> TARGET. In-place: SORTED -> TARGET (== INPUT).
TARGET_UNCOMPRESSED="${TARGET%.gz}"
mv "$SORTED" "$TARGET_UNCOMPRESSED"
rm -f "$TMP"
# Remove any pre-existing .gz at the target (in in-place mode the original
# plain-gzip file is still present; bgzip refuses to overwrite without -f).
rm -f "$TARGET"
bgzip -@ 4 "$TARGET_UNCOMPRESSED"
echo "[$(date '+%F %T')] bgzipped: $TARGET ($(stat -c %s "$TARGET") bytes)"

rm -f "${TARGET}.tbi"
tabix -p "$FMT" "$TARGET"
echo "[$(date '+%F %T')] indexed: ${TARGET}.tbi ($(stat -c %s "${TARGET}.tbi") bytes)"

# Sanity check: pull the first contig's first 100 kb and confirm tabix returns rows.
FIRST_CONTIG=$(zcat "$TARGET" | awk '$1!~/^#/ {print $1; exit}')
if [[ -n "$FIRST_CONTIG" ]]; then
    N=$(tabix "$TARGET" "${FIRST_CONTIG}:1-100000" | wc -l)
    echo "[$(date '+%F %T')] sanity: ${FIRST_CONTIG}:1-100000 returns $N row(s)"
fi

if [[ "$TARGET" == "$INPUT" ]]; then
    echo "[$(date '+%F %T')] DONE — track ready for igv-reports. Original preserved at $BACKUP"
else
    echo "[$(date '+%F %T')] DONE — sibling track ready at $TARGET. Original $INPUT untouched."
fi
