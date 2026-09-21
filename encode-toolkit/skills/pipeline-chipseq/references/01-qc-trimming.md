# Stage 1: QC and Trimming

## Input
- Raw FASTQ files (single-end or paired-end)
- Adapter sequences (auto-detected by Trim Galore or specify: Illumina TruSeq)

## Tools
- **FastQC 0.12.1** (image version): Per-base quality, adapter content, duplication rates, GC content
- **Trim Galore 0.6.10** (image version, wraps Cutadapt 4.6): Adapter trimming + quality filtering

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| Quality cutoff | 20 | Phred score minimum |
| Min length | 36 | Discard reads shorter than this after trimming |
| Adapter | auto-detect | TruSeq for most ENCODE libraries |
| Stringency | 1 | Overlap with adapter sequence required |
| Error rate | 0.1 | Maximum allowed error rate in adapter detection |

The workflow does not expose these as parameters; they are fixed in `main.nf`.

## Commands

The workflow runs the equivalent of:

```bash
# Raw QC
fastqc -t 4 -o qc_raw/ sample_R1.fastq.gz sample_R2.fastq.gz

# Paired-end trimming
trim_galore --paired --quality 20 --length 36 --fastqc \
  --cores 4 -o trimmed/ sample_R1.fastq.gz sample_R2.fastq.gz

# Single-end trimming (--single_end)
trim_galore --quality 20 --length 36 --fastqc \
  --cores 4 -o trimmed/ sample.fastq.gz
```

## Expected Output
- `trimmed/*_trimming_report.txt` -- trimming statistics (reads processed, trimmed, removed)
- `trimmed/*_val_1.fq.gz`, `trimmed/*_val_2.fq.gz` -- trimmed paired-end reads
- `trimmed/*_trimmed.fq.gz` -- trimmed single-end reads
- `fastqc/*.html`, `fastqc/*.zip` -- FastQC reports for raw reads (FASTQC process) and for
  the trimmed reads (`trim_galore --fastqc`). Both land in the same `fastqc/` directory;
  the trimmed-read reports carry `_val_1`/`_val_2`/`_trimmed` in the file name.

All of these feed the MultiQC report.

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| Per-base quality | >Q20 across all positions after trimming | Check sequencing run quality |
| Adapter content | <5% after trimming | Verify trimming parameters |
| GC content | Unimodal, matching expected genome GC | Check for contamination |
| Sequence duplication | <50% at this stage | May indicate low complexity library |
| Read count | Record for downstream normalization | No hard threshold at this stage |

## Troubleshooting

- **High adapter content after trimming**: the workflow does not expose `--stringency` or
  `--adapter`; run Trim Galore manually with those flags if auto-detection fails
- **Bimodal GC distribution**: Indicates possible contamination; run FastQ Screen to identify
  organism of origin
- **Very short reads after trimming**: Library insert size may be too short; consider
  investigating library preparation (the 36 bp minimum is fixed in `main.nf`)
