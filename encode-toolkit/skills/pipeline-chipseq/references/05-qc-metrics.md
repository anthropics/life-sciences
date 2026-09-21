# Stage 5: Signal Tracks and QC Report

Stage 5 of the workflow does three things: it converts the MACS2 bedGraphs into bigWig
signal tracks, it computes FRiP for every treatment sample, and it runs MultiQC over the
logs collected in earlier stages. The sections marked "manual step" below are **not run by
the workflow**.

## Signal Track Generation (workflow)

Signal tracks provide normalized coverage for genome browser visualization. The workflow
runs the equivalent of the following per sample, using the `_treat_pileup.bdg` and
`_control_lambda.bdg` files MACS2 wrote with `-B`. MACS2 produces `_control_lambda.bdg`
even without `-c`, so the tracks exist with or without `--control`; without a control the
"fold change" track is fold enrichment over the local lambda model.

```bash
# Fold enrichment over the MACS2 control/lambda track
macs2 bdgcmp -t sample_treat_pileup.bdg -c sample_control_lambda.bdg \
  -o fc.bdg -m FE
sort -k1,1 -k2,2n fc.bdg > fc.sorted.bdg
bedGraphToBigWig fc.sorted.bdg GRCh38.chrom.sizes sample.fc.bw

# Signal p-value track (statistical significance)
macs2 bdgcmp -t sample_treat_pileup.bdg -c sample_control_lambda.bdg \
  -o pval.bdg -m ppois
sort -k1,1 -k2,2n pval.bdg > pval.sorted.bdg
bedGraphToBigWig pval.sorted.bdg GRCh38.chrom.sizes sample.pval.bw
```

Published as `signal/<sample>.fc.bw` and `signal/<sample>.pval.bw`. `--chrom_sizes` is
required for this stage, which is why the workflow stops at startup when it is missing.

## FRiP (workflow)

The workflow computes the fraction of reads in peaks for every treatment sample, in both
narrow and broad mode, and publishes `qc/<sample>.frip_mqc.tsv`. It runs the equivalent of:

```bash
TOTAL_READS=$(samtools view -c results/filtered/sample.final.bam)
READS_IN_PEAKS=$(bedtools intersect -u -a results/filtered/sample.final.bam \
  -b results/peaks/narrow/sample_peaks.narrowPeak | samtools view -c -)
awk -v a=$READS_IN_PEAKS -v b=$TOTAL_READS 'BEGIN{printf "%.4f\n", a/b}'
```

The division uses `awk`. The published file carries a MultiQC header and one row per peak
file, so MultiQC renders it as the section "Fraction of reads in peaks":

| Column | Meaning |
|--------|---------|
| `Peak set` | Name of the peak file the row was computed from |
| `FRiP` | `reads_in_peaks / total_reads`, printed to four decimals |
| `reads_in_peaks` | Alignments of the final BAM overlapping a peak (`bedtools intersect -u`) |
| `total_reads` | All alignments of the final BAM (`samtools view -c`) |

FRiP here is a fraction, not a percentage: the ENCODE >=1% standard is >=0.01 in this
file. Control libraries are not peak-called, so no `CONTROL_<name>.frip_mqc.tsv` is
written.

## MultiQC Aggregated Report (workflow)

```bash
multiqc . -o . -f
```

Published as `qc/multiqc/multiqc_report.html` with `qc/multiqc/multiqc_data/`. The workflow
feeds it exactly these inputs:

- FastQC reports for raw reads
- FastQC reports for trimmed reads (from `trim_galore --fastqc`)
- Trim Galore trimming reports
- `samtools flagstat` from the alignment step
- Picard MarkDuplicates metrics
- `samtools flagstat` from the blacklist-filtered BAM
- `qc/<sample>.frip_mqc.tsv` from the FRiP step, shown as "Fraction of reads in peaks"

MACS2 output is not passed to MultiQC.

## Comprehensive QC Metrics Table

| Metric | Tool | Threshold | Computed by the workflow? |
|--------|------|-----------|---------------------------|
| Total reads | FastQC, samtools flagstat | Record | yes |
| Mapped reads | samtools flagstat | >=20M TF / >=45M histone | yes |
| Mapping rate | samtools flagstat | >80% | yes |
| Duplication rate | Picard | <30% | yes |
| IDR peaks | IDR | >20,000 (TF) | yes (narrow runs, one file per replicate pair) |
| NRF | manual (see 03-filtering.md) | >=0.8 | no |
| PBC1 | manual (see 03-filtering.md) | >=0.8 | no |
| PBC2 | manual (see 03-filtering.md) | >=3 | no |
| NSC | phantompeakqualtools | >1.05 | no |
| RSC | phantompeakqualtools | >0.8 | no |
| FRiP | bedtools + samtools (`qc/<sample>.frip_mqc.tsv`) | >=0.01 (1%) | yes |
| Mitochondrial fraction | `samtools idxstats` | <5% | no |

## Manual step: Strand Cross-Correlation (phantompeakqualtools)

Computes NSC and RSC, which measure ChIP enrichment quality independent of peak calling.
Reference: Kharchenko et al. 2008; Landt et al. 2012.

**Not run by this workflow, and phantompeakqualtools is not in the pipeline image.** Run it
in the conda environment `bioinformatics-installer/environments/chipseq-env.yml`, which
includes `phantompeakqualtools=1.2.2`.

```bash
conda env create -f chipseq-env.yml && conda activate encode-chipseq
run_spp.R -c=results/filtered/sample.final.bam -savp=cc_plot.pdf -out=cc_scores.txt -tmpdir=tmp/
# Output columns: filename, numReads, estFragLen, corr_estFragLen,
#                 phantomPeak, corr_phantomPeak, argmin_corr, min_corr,
#                 NSC, RSC, QualityTag
```

## Manual step: Fingerprint Plot (deeptools)

**Not run by this workflow**, though deeptools 3.5.5 is installed in the image:

```bash
plotFingerprint -b results/filtered/chip.final.bam results/filtered/CONTROL_input.final.bam \
  --labels ChIP Input \
  --plotFile fingerprint.pdf \
  --outRawCounts fingerprint_counts.txt
```

The fingerprint plot visually shows enrichment: a good ChIP sample curves away
from the diagonal (uniform coverage), while input stays close to the diagonal.
