# Stage 1: QC and Trimming

## Input
- Raw paired-end FASTQ files (this workflow is paired-end only)
- Adapter sequences: Nextera transposase adapters (not TruSeq)

## Tools
- **FastQC 0.12.1** (image version): Per-base quality, adapter content, duplication rates, insert size
- **Trim Galore 0.6.10** (image version, wraps Cutadapt 4.6): Adapter trimming + quality filtering

## Key Difference from ChIP-seq
ATAC-seq uses **Nextera** transposase adapters, not Illumina TruSeq. The workflow passes
`--nextera` explicitly rather than relying on auto-detection.

## Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| Quality cutoff | 20 | Phred score minimum |
| Min length | 20 | Shorter than ChIP-seq due to NFR fragments |
| Adapter | Nextera (`--nextera`) | Tn5 transposase adapters |
| Stringency | 1 | Overlap with adapter sequence required |

The workflow does not expose these as parameters; they are fixed in `main.nf`.

## Commands

The workflow runs the equivalent of:

```bash
# Raw QC
fastqc -t 4 -o qc_raw/ sample_R1.fastq.gz sample_R2.fastq.gz

# Paired-end trimming with Nextera adapters
trim_galore --paired --nextera --quality 20 --length 20 --fastqc \
  --cores 4 -o trimmed/ sample_R1.fastq.gz sample_R2.fastq.gz
```

## Expected Output
- `trimmed/*_trimming_report.txt` -- trimming statistics
- `trimmed/*_val_1.fq.gz`, `trimmed/*_val_2.fq.gz` -- trimmed paired-end reads
- `fastqc/*.html`, `fastqc/*.zip` -- FastQC reports for raw reads (FASTQC process) and for
  the trimmed reads (`trim_galore --fastqc`). Both land in the same `fastqc/` directory;
  the trimmed-read reports carry `_val_1`/`_val_2` in the file name.

All of these feed the MultiQC report.

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| Per-base quality | >Q20 after trimming | Check sequencing run quality |
| Adapter content | <5% after trimming | Verify Nextera adapter detection |
| GC content | Unimodal, matching genome | Check for contamination |
| Insert size | Nucleosomal ladder visible | Inspect library prep |
| Read count | >=50M total recommended | May need more sequencing |

## Troubleshooting
- **High adapter content**: ATAC-seq libraries with short inserts have more adapter
  read-through. This is normal for NFR fragments; the fixed `--length 20` keeps them.
- **No nucleosomal pattern in insert sizes**: May indicate failed transposition or
  over-transposition. Check Tn5:cell ratio.
