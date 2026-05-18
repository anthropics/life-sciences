#!/usr/bin/env bash
# examples/portable/single_sample.sh — generic single-sample build.
#
# Builds one HTML for a tumor/normal pair at a handful of SNV sites.
# Assumes:
#   - `pip install igv-reports` has put `create_report` on PATH.
#   - You have your own hg38 FASTA (with .fai sibling) and BAMs.
#   - You have your own gencode + CpG-islands track files (or skip them
#     with --no-default-tracks alone).
#
# Set these to match your environment before running.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORKDIR="${WORKDIR:-${PWD}/igv_reports_demo}"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

# --- inputs (edit these) ---
FASTA="${FASTA:-${HOME}/data/hg38/hg38.fa}"          # must have ${FASTA}.fai
TUMOR_BAM="${TUMOR_BAM:-${HOME}/data/tumor.bam}"
NORMAL_BAM="${NORMAL_BAM:-${HOME}/data/normal.bam}"
GENCODE_GFF="${GENCODE_GFF:-${HOME}/data/hg38/gencode.v47.annotation.gff3.gz}"  # bgzip+tabix
CPG_ISLANDS="${CPG_ISLANDS:-${HOME}/data/hg38/hg38_CpGIslands.bed}"

# --- sites BED (4 cols: chrom, start, end, name) ---
cat > sites.hg38.bed <<'EOF'
#chrom	start	end	name
chr2	25246499	25246500	DNMT3A_R882
chr7	148884000	148884001	EZH2_Y646
EOF

python "${REPO_ROOT}/scripts/build_igvreports.py" \
    --genome hg38 \
    --sites sites.hg38.bed \
    --bam "${TUMOR_BAM}" "${NORMAL_BAM}" \
    --fasta "${FASTA}" \
    --no-default-tracks \
    --extra-track "${GENCODE_GFF}" \
    --extra-track "${CPG_ISLANDS}" \
    --info-columns name \
    --output report.hg38.html \
    --no-apptainer

echo "Done. Open ${WORKDIR}/report.hg38.html in a browser."
