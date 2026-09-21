# Stage 5: QC Metrics

## Tools
- **RSeQC 5.0.3**: RNA-seq quality control suite (Wang et al. 2012, ~3,500 citations)
- **MultiQC 1.21**: Aggregated QC report generation

The workflow runs four RSeQC modules — `infer_experiment.py`, `read_distribution.py`,
`geneBody_coverage.py` and, for paired-end runs, `inner_distance.py` — all against the
BED12 file given by `--rseqc_bed`, and publishes their output to `qc/rseqc/`. R is not
installed in the image, so the modules that end by calling Rscript write their `.txt` and
`.r` files but no rendered plot.

## RSeQC Modules

### infer_experiment.py (Strandedness Check)

```bash
infer_experiment.py -r hg38_RefSeq.bed -i sample.Aligned.sortedByCoord.out.bam \
  > sample.infer_experiment.txt
```

Determines library strandedness by sampling read orientation relative to annotated
transcripts. This is a post-hoc check: it does not feed back into quantification. If it
disagrees with the `--strandedness` value the run used, rerun the pipeline with the
correct value.

| Output Pattern | Interpretation | `--strandedness` to use |
|----------------|---------------|-------------|
| "1++,1--,2+-,2-+" > 90% | Forward stranded | `forward` |
| "1+-,1-+,2++,2--" > 90% | Reverse stranded (dUTP) | `reverse` |
| ~50/50 split | Unstranded | `none` |

**ENCODE standard**: Expect >90% reverse-stranded reads for dUTP libraries.

### read_distribution.py (Mapping Distribution)

```bash
read_distribution.py -r hg38_RefSeq.bed -i sample.Aligned.sortedByCoord.out.bam \
  > sample.read_distribution.txt
```

Reports fraction of reads mapping to CDS exons, 5' UTR, 3' UTR, introns, and
intergenic regions.

| Region | Expected (mRNA-seq) | Concern Threshold |
|--------|--------------------|--------------------|
| CDS exons | 40-60% | <30% suggests degradation or DNA contamination |
| 5' UTR | 5-10% | <2% suggests 5' degradation |
| 3' UTR | 15-25% | >40% suggests 3' bias (degraded RNA) |
| Introns | 10-25% | >40% suggests DNA contamination or pre-mRNA |
| Intergenic | <5% | >10% suggests DNA contamination |

### geneBody_coverage.py (Gene Body Coverage)

```bash
geneBody_coverage.py -r hg38_RefSeq.bed \
  -i sample.Aligned.sortedByCoord.out.bam -o sample.geneBody_coverage
```

The workflow passes the same `--rseqc_bed` file used by the other modules. A
housekeeping-gene BED (for example `hg38_HouseKeeping.bed`) runs faster and is the usual
choice when calling this module by hand; both are valid inputs.

Reports normalized coverage across gene bodies (5' to 3') in
`<sample>.geneBody_coverage.geneBodyCoverage.txt`. Uniform coverage indicates intact RNA;
strong 3' bias indicates degradation.

| Pattern | Interpretation |
|---------|---------------|
| Uniform (5'/3' ratio 0.7-1.3) | Good RNA quality |
| 3' bias (5'/3' ratio <0.5) | RNA degradation |
| 5' bias (5'/3' ratio >2.0) | Possible oligo-dT priming bias |

### inner_distance.py (Insert Size Distribution)

```bash
inner_distance.py -r hg38_RefSeq.bed \
  -i sample.Aligned.sortedByCoord.out.bam -o sample.inner_distance
```

Paired-end only; skipped when the run uses `--single_end`. Reports the inner distance
between mates across a set of files matching `<sample>.inner_distance.*`. Negative values
indicate overlapping reads (common for short inserts). The peak should match the expected
library insert size (typically 150-300 bp).

### STAR Log Metrics

The `star/<sample>.Log.final.out` from STAR provides critical metrics:

| Metric | Threshold | Notes |
|--------|-----------|-------|
| Uniquely mapped reads % | >=70% | Primary quality indicator |
| Multi-mapped reads % | <10% | High suggests repetitive contamination |
| Unmapped: too short % | <10% | High suggests over-trimming |
| % of reads mapped to multiple loci | <10% | Expected for gene families |
| % of chimeric reads | <1% | Zero unless STAR is run with chimeric detection |
| Number of splices: Total | Millions expected | Low count suggests annotation mismatch |

## Manual Checks (not run by this workflow)

The metrics below are not produced by the pipeline. Run them yourself against the
published BAM and RSEM output when you need them, and label them as separate steps when
reporting.

### rRNA Rate Assessment

```bash
# Count mapped reads overlapping rRNA loci. -L takes a BED file of regions; without it
# samtools would read the file name as a region string and fail.
samtools view -c -F 4 -L rRNA_intervals.bed star/sample.Aligned.sortedByCoord.out.bam
```

Divide by the total mapped count from `Log.final.out` to get the rate.

| rRNA Rate | Interpretation |
|-----------|---------------|
| <5% | Excellent rRNA depletion |
| 5-10% | Acceptable |
| 10-30% | Suboptimal; reduced effective depth |
| >30% | Failed rRNA depletion; consider re-prep |

### Saturation Analysis

```bash
# RSeQC RPKM saturation (installed in the image, but never invoked by the workflow;
# it writes an .r script that needs R to render)
RPKM_saturation.py -r hg38_RefSeq.bed \
  -i star/sample.Aligned.sortedByCoord.out.bam -o sample_saturation
```

Subsamples reads at increasing fractions (5%, 10%, ..., 100%) and measures gene
detection. A plateau indicates sufficient sequencing depth. If the curve is still
rising at 100%, more sequencing is recommended.

### Detected Genes

```bash
# Genes with TPM > 1; the TPM column of RSEM genes.results is column 6
awk 'NR > 1 && $6 > 1 {n++} END {print "Genes with TPM>1:", n}' rsem/sample.genes.results
```

Expect >12,000 for a human sample at ENCODE depth.

### Library Duplication Rate

Picard `MarkDuplicates` and `CollectRnaSeqMetrics` are not in the container image. The
FastQC report does carry a sequence-level duplication estimate, which is not the same
quantity; run Picard separately if you need the alignment-based rate.

## MultiQC Aggregation

```bash
multiqc . -o . -f
```

The workflow feeds MultiQC exactly these inputs: FastQC reports for the raw reads, the
Trim Galore trimming reports, the FastQC reports for the trimmed reads, the STAR
`Log.final.out`, the RSEM `.stat/` directories, and the RSeQC `infer_experiment` and
`read_distribution` outputs. `geneBody_coverage`, `inner_distance` and the Kallisto
output are published but are not part of the report — read those files directly. The
report is written to `qc/multiqc/multiqc_report.html` with the parsed values in
`qc/multiqc/multiqc_data/`.

## Expected Output
- `qc/rseqc/<sample>.infer_experiment.txt` -- strandedness inference
- `qc/rseqc/<sample>.read_distribution.txt` -- mapping distribution by genomic feature
- `qc/rseqc/<sample>.geneBody_coverage.geneBodyCoverage.txt` and `.r` -- gene body coverage
  (RSeQC also tries to draw `.curves.pdf`, which needs R and is therefore absent)
- `qc/rseqc/<sample>.inner_distance.*` -- insert size, paired-end runs only (RSeQC also
  tries to draw a PDF, which needs R)
- `qc/multiqc/multiqc_report.html` -- aggregated QC report
