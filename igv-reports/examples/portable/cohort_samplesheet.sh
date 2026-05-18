#!/usr/bin/env bash
# examples/portable/cohort_samplesheet.sh — generic cohort build.
#
# Builds one HTML per row of a TSV samplesheet, plus an index.html linking
# them all. Demonstrates the samplesheet format and the most common flags.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORKDIR="${WORKDIR:-${PWD}/igv_reports_cohort_demo}"
mkdir -p "$WORKDIR" && cd "$WORKDIR"

# --- inputs (edit these) ---
FASTA="${FASTA:-${HOME}/data/hg38/hg38.fa}"
GENCODE_GFF="${GENCODE_GFF:-${HOME}/data/hg38/gencode.v47.annotation.gff3.gz}"

# --- samplesheet (one row per sample) ---
# Required columns: sample, sites_bed
# Optional columns: bam_tumor, bam_normal, vcf, extra_tracks (comma-separated)
cat > cohort.tsv <<EOF
sample	bam_tumor	bam_normal	vcf	sites_bed
p001	${HOME}/data/p001/tumor.bam	${HOME}/data/p001/normal.bam		sites.hg38.bed
p002	${HOME}/data/p002/tumor.bam	${HOME}/data/p002/normal.bam		sites.hg38.bed
EOF

# --- shared sites for both patients ---
cat > sites.hg38.bed <<'EOF'
#chrom	start	end	name
chr2	25246499	25246500	DNMT3A_R882
chr7	148884000	148884001	EZH2_Y646
EOF

python "${REPO_ROOT}/scripts/build_igvreports.py" \
    --genome hg38 \
    --samplesheet cohort.tsv \
    --fasta "${FASTA}" \
    --no-default-tracks \
    --extra-track "${GENCODE_GFF}" \
    --output-dir reports \
    --no-apptainer

echo "Done. Open ${WORKDIR}/reports/index.html in a browser."
echo "Cohort verifier ran automatically; see reports/cohort_verify.summary.md."
