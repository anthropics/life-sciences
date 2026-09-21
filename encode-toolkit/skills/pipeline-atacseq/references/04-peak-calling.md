# Stage 4: Peak Calling and IDR

## Tools
- **MACS2 2.2.9.1** (image version): Peak caller (Zhang et al. 2008)
- **IDR 2.0.4.2** (image version): Irreproducible Discovery Rate (Li et al. 2011)
- **bedtools 2.31.0 + samtools 1.19** (image versions): FRiP counting after peak calling

## MACS2 Parameters for ATAC-seq

ATAC-seq peak calling differs from ChIP-seq in several key ways:
1. **No control/input file** -- peaks are called against the local background model
2. **Use the NFR BAM only** -- nucleosome-free fragments for accessibility peaks
3. **BAMPE mode** -- fragment coordinates come from the read pairs

| Parameter | Value | Notes |
|-----------|-------|-------|
| `--format` | BAMPE | Use actual fragment sizes |
| `--gsize` | hs (2.7e9) | `mm` when `--genome mm10` |
| `--nomodel` | yes | Do not build a shifting model |
| `--qvalue` | 0.05 | FDR threshold |
| `--keep-dup` | all | Duplicates were already removed in Stage 3 |
| `--call-summits` | yes | Identify sub-peak summits |
| `-B` | yes | Generate bedGraphs |

**`--shift` and `--extsize` are not used, and would have no effect here.** In `-f BAMPE`
mode MACS2 takes fragment coordinates from the read pairs: it forces `nomodel = True` and
sets `shift = 0` internally, without warning. Shift/extension values such as
`--shift -100 --extsize 200` or `--shift -75 --extsize 150` apply only when calling peaks
on BED or single-end input, which this workflow does not do.

## Commands

The workflow runs the equivalent of:

```bash
# Peak calling on NFR fragments
macs2 callpeak -t sample.nfr.bam \
  -f BAMPE -g hs -n sample \
  --nomodel --keep-dup all --call-summits \
  --qvalue 0.05 -B

# IDR, once per pair of replicates, the two names sorted alphabetically
idr --samples rep1_peaks.narrowPeak rep2_peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file rep1_vs_rep2.idr_peaks.txt \
  --plot \
  --idr-threshold 0.05
```

Calling peaks on the blacklist-filtered all-fragment BAM (`filtered/<sample>.final.bam`)
instead of the NFR BAM is a reasonable manual variant, but the workflow always uses NFR.

## What the workflow publishes

- `peaks/narrow/<sample>_peaks.narrowPeak`
- `peaks/narrow/<sample>_peaks.xls`
- `peaks/narrow/<sample>_treat_pileup.bdg` and `<sample>_control_lambda.bdg`
- `peaks/idr/<sampleA>_vs_<sampleB>.idr_peaks.txt` and, when IDR emits it,
  `<sampleA>_vs_<sampleB>.idr_peaks.txt.png`, one pair of files per replicate pair
- `qc/<sample>.frip_mqc.tsv` (see `references/05-qc-metrics.md`)

`<sample>_summits.bed` is produced by `--call-summits` and published to `peaks/narrow/`
alongside the peak files. Only narrow peaks are called; there is no broad mode.

## What the IDR step does and does not do

All samples matched by `--reads` are treated as replicates of one experiment. The workflow
sorts them by sample name and runs `idr` once for every pair. Consequently:

- Two samples give one comparison, three give three, four give six. No replicate is
  dropped.
- With a single peak file there is no pair, so IDR is skipped silently and `peaks/idr/` is
  not created.
- There are **no** pooled-replicate calls, **no** pseudoreplicates, and **no** optimal /
  conservative peak sets.
- Rescue ratio and self-consistency ratio are **not computed**; they require pooled and
  pseudoreplicated peak calls that this workflow does not produce.

`--skip_idr` turns the step off entirely.

## IDR Interpretation

| Metric | Expected (Good) | Concern | Computed here? |
|--------|-----------------|---------|----------------|
| IDR peaks (0.05 threshold) | 50,000-150,000 | <30,000 suggests poor signal | yes |
| Rescue ratio | <2 | >2 suggests replicate discordance | no |
| Self-consistency ratio | <2 | >2 suggests noisy data | no |

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| FRiP (`qc/<sample>.frip_mqc.tsv`) | >=0.3 (ATAC-seq standard) | Poor accessibility signal |
| Peak count | >50,000 (IDR filtered) | Low enrichment |
| Peak width distribution | Median 200-500 bp | Check if calling mode correct |
| Peaks at TSS (manual) | Enrichment visible | Fundamental ATAC-seq signal |

## Notes

- ATAC-seq FRiP is typically much higher than ChIP-seq (0.3-0.6 vs 0.01-0.1)
  because open chromatin is a large fraction of the genome. The workflow computes it after
  peak calling, counting the peaks called here against the blacklist-filtered all-fragment
  BAM; see `references/05-qc-metrics.md`.
- The workflow always calls peaks on NFR fragments for accessibility analysis.
- For nucleosome positioning, use `filtered/nfr/<sample>.mononuc.bam` separately.
- IDR is standard for ATAC-seq with biological replicates.
