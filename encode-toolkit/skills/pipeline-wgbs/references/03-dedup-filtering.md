# Deduplication and Filtering for WGBS

PCR duplicates inflate coverage estimates and bias methylation calls.
Deduplication is essential for WGBS but must NOT be used for RRBS.

## Bismark Deduplication

Bismark provides its own deduplication tool optimized for bisulfite data, and this is what
the workflow runs (skipped with `--skip_dedup`):

```bash
deduplicate_bismark \
    --bam \
    --paired \
    bismark_out/sample_pe.bam
```

Bismark deduplication identifies duplicates by their alignment positions on
both strands, accounting for bisulfite conversion. This is preferred over
Picard for bisulfite data; Picard is not installed in the image.

The deduplicated BAM is then sorted and indexed with samtools and published as
`bismark/alignments/<sample>.sorted.bam` (+ `.bai`). The deduplication report goes to
`bismark/dedup_reports/`.

### Expected Duplication Rates

| Library Quality | Duplication Rate | Action |
|----------------|------------------|--------|
| Good | <20% | Proceed |
| Acceptable | 20-40% | Proceed with caution |
| Poor | 40-60% | Consider resequencing |
| Very poor | >60% | Library failed -- redo |

## Manual Alternative: Picard MarkDuplicates (not run by this workflow)

Picard is neither in the container image nor used by the workflow. If you align with
bwa-meth by hand (`02-bismark-alignment.md`), which has no built-in dedup, install Picard
separately and run:

```bash
picard MarkDuplicates \
    INPUT=sample_bwameth_sorted.bam \
    OUTPUT=sample_dedup.bam \
    METRICS_FILE=sample_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true
```

## Optional Manual BAM Filtering (not run by this workflow)

The workflow hands MethylDackel the deduplicated, sorted BAM with no MAPQ or flag
filtering. If you want a filtered BAM, produce it yourself from the published
`<sample>.sorted.bam` and point MethylDackel at that instead:

```bash
samtools view -b -h \
    -q 10 \
    -F 1804 \
    -f 2 \
    bismark/alignments/sample.sorted.bam \
    | samtools sort -@ 4 -o sample_filtered.bam

samtools index sample_filtered.bam
```

### Filter Flag Explanation

| Flag | Binary | Meaning |
|------|--------|---------|
| `-q 10` | MAPQ >= 10 | Minimum mapping quality |
| `-F 4` | 0x4 | Remove unmapped reads |
| `-F 8` | 0x8 | Remove reads whose mate is unmapped |
| `-F 256` | 0x100 | Remove secondary alignments |
| `-F 512` | 0x200 | Remove reads failing QC |
| `-F 1024` | 0x400 | Remove PCR duplicates |
| `-f 2` | 0x2 | Keep only properly paired |

Combined: `-F 1804` = 4 + 8 + 256 + 512 + 1024, removing unmapped, mate-unmapped,
secondary, QC-fail and duplicate reads.

## RRBS: Skip Deduplication

For RRBS libraries, MspI digestion creates identical fragment starts at cut
sites. These are NOT PCR duplicates and must be retained. The trimming in this workflow is
WGBS-specific (`01-qc-trimming.md`), so RRBS is better processed by hand; if you run it
through the workflow anyway, at least set `--skip_dedup true`.

```bash
# Do NOT run deduplication for RRBS.
# To filter by hand instead, drop the duplicate bit from the mask (780 = 4+8+256+512):
samtools view -b -h -q 10 -F 780 -f 2 sample.bam \
    | samtools sort -@ 4 -o sample_filtered.bam
```

## Coverage Statistics

The workflow computes its coverage statistics from the CpG bedGraph, not from the BAM;
see `05-qc-metrics.md`. The commands below are alternatives you can run by hand on the
published BAM. `bedtools` is not in the container image — it ships with the conda
environment instead.

```bash
samtools depth -a bismark/alignments/sample.sorted.bam \
    | awk '{sum+=$3; n++} END {print "Mean coverage:", sum/n}'

samtools flagstat bismark/alignments/sample.sorted.bam > sample_flagstat.txt

# CpG-specific coverage: reads overlapping each CpG (last column), as a histogram.
# The BED goes in -a so the output is text; with a BAM in -a, bedtools writes BAM.
bedtools coverage \
    -a /ref/CpG_sites.bed \
    -b bismark/alignments/sample.sorted.bam \
    -counts \
    | awk '{print $NF}' | sort -n | uniq -c
```
