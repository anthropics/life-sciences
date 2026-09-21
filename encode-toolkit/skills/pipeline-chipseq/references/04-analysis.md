# Stage 4: Peak Calling and IDR

## Tools
- **MACS2 2.2.9.1** (image version): Model-based Analysis of ChIP-Seq (Zhang et al. 2008, ~7,000 citations)
- **IDR 2.0.4.2** (image version): Irreproducible Discovery Rate (Li et al. 2011, ~1,500 citations)
- **bedtools 2.31.0 + samtools 1.19** (image versions): FRiP counting after peak calling

## MACS2 Parameters

The workflow builds one MACS2 command per ChIP sample. `--peak_type` selects the narrow or
broad column for the whole run.

| Parameter | Narrow (TF) | Broad (Histone) | Notes |
|-----------|-------------|-----------------|-------|
| `--format` | BAMPE | BAMPE | `BAM` when `--single_end` is set |
| `--gsize` | hs (2.7e9) | hs (2.7e9) | `mm` when `--genome mm10` |
| `--qvalue` | 0.05 | 0.05 | FDR threshold |
| `--broad` | no | yes | Broad peak mode for repressive marks |
| `--broad-cutoff` | n/a | 0.1 | Linking threshold for broad peaks |
| `--nomodel` | yes | yes | Fragment sizes come from the BAM |
| `--keep-dup` | all | all | Duplicates were already removed in Stage 3 |
| `--call-summits` | yes | no | Subpeak summit positions |
| `-B` | yes | yes | bedGraphs for the Stage 5 signal tracks |
| `-c` | only with `--control` | only with `--control` | All control BAMs are passed together; MACS2 pools them |

## Narrow vs Broad Mark Decision

| Peak Type | Targets | Rationale |
|-----------|---------|-----------|
| **Narrow** | H3K4me3, H3K4me1, H3K27ac, H3K9ac, all TFs, CTCF | Punctate binding pattern |
| **Broad** | H3K27me3, H3K36me3, H3K9me3, H3K79me2 | Diffuse domain spreading |

## Commands

The workflow runs the equivalent of:

```bash
# Narrow peaks (TF and active histone marks)
macs2 callpeak -t treatment.bam -c control.bam \
  -f BAMPE -g hs -n sample \
  --qvalue 0.05 --nomodel --keep-dup all --call-summits -B

# Broad peaks (repressive histone marks)
macs2 callpeak -t treatment.bam -c control.bam \
  -f BAMPE -g hs -n sample \
  --qvalue 0.05 --nomodel --keep-dup all --broad --broad-cutoff 0.1 -B

# IDR, once per pair of replicates, the two names sorted alphabetically (narrow runs only)
idr --samples rep1_peaks.narrowPeak rep2_peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file rep1_vs_rep2.idr_peaks.txt \
  --plot \
  --idr-threshold 0.05
```

## FRiP (workflow)

After peak calling the workflow computes the fraction of reads in peaks for every
treatment sample, in both narrow and broad mode, and publishes
`qc/<sample>.frip_mqc.tsv`. The command and the output columns are in
`references/05-qc-metrics.md`.

## What the workflow publishes

- `peaks/<narrow|broad>/<sample>_peaks.narrowPeak` (or `_peaks.broadPeak` and
  `_peaks.gappedPeak` in broad mode)
- `peaks/<narrow|broad>/<sample>_peaks.xls`
- `peaks/<narrow|broad>/<sample>_treat_pileup.bdg` and `<sample>_control_lambda.bdg`
  (consumed by Stage 5)
- `peaks/idr/<sampleA>_vs_<sampleB>.idr_peaks.txt` and, when IDR emits it,
  `<sampleA>_vs_<sampleB>.idr_peaks.txt.png`, one pair of files per replicate pair
- `qc/<sample>.frip_mqc.tsv` (see `references/05-qc-metrics.md`)

`<sample>_summits.bed` is produced by `--call-summits` and published to `peaks/narrow/`;
broad runs do not produce it.

## What the IDR step does and does not do

All samples matched by `--reads` are treated as replicates of one experiment. The workflow
sorts them by sample name and runs `idr` once for every pair. Consequently:

- Only narrow runs get IDR (`--peak_type broad` produces no `peaks/idr/`).
- Two samples give one comparison, three give three, four give six. No replicate is
  dropped.
- With a single peak file there is no pair, so IDR is skipped silently.
- There are **no** pooled-replicate calls, **no** pseudoreplicates, and **no** optimal /
  conservative peak sets.
- Rescue ratio and self-consistency ratio are **not computed**. Computing them requires
  pooled and pseudoreplicated peak calls that this workflow does not produce.

## Manual follow-ups (not run by this workflow)

```bash
# Self-pseudoreplicate IDR: split one replicate's reads in half, call peaks on each half,
# then compare. Requires re-running MACS2 outside the workflow.
idr --samples pr1_peaks.narrowPeak pr2_peaks.narrowPeak \
  --input-file-type narrowPeak \
  --output-file self_idr_peaks.txt \
  --plot
```

Pooled peak calls across replicates are likewise outside this workflow.

## IDR Interpretation

| Metric | Expected (Good) | Concern | Computed here? |
|--------|-----------------|---------|----------------|
| IDR peaks (0.05 threshold) | 50,000-200,000 (TF) | <20,000 suggests poor enrichment | yes |
| IDR / individual rep ratio | 0.3-0.7 | <0.3 very stringent; >0.7 very lenient | manual (count lines) |
| Rescue ratio | <2 | >2 suggests replicate discordance | no |
| Self-consistency ratio | <2 | >2 suggests noisy data | no |

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| FRiP (`qc/<sample>.frip_mqc.tsv`) | >=0.01 (1%) | Poor enrichment; check antibody |
| IDR peaks | >20,000 (TF) | Low enrichment or poor replicates |
| Peak count ratio | IDR = 30-70% of rep peaks | Replicate consistency issue |

## Notes

- Peaks are called on individual replicates; the IDR step then compares every pair of them.
- For broad marks, IDR is not standard practice, which matches the workflow: broad runs
  skip it. Use pooled replicate peak calls instead, produced manually.
- The `--call-summits` flag identifies subpeak summits within broader peak regions,
  useful for motif analysis downstream.
- FRiP is calculated by the workflow as: alignments of `filtered/<sample>.final.bam` that
  overlap a peak / all alignments of that BAM. Control libraries are not peak-called, so
  they get no FRiP row.
