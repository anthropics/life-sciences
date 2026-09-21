# Stage 3: Filtering and Deduplication

## Tools
- **Samtools 1.19** (image version): Flag and MAPQ filtering
- **Picard MarkDuplicates 3.1.1** (image version, Java 17): Remove PCR duplicates
- **bedtools 2.31.0** (image version): Blacklist region filtering

## Order of operations in the workflow

```
aligned.bam -> samtools flag/MAPQ filter -> Picard MarkDuplicates (REMOVE_DUPLICATES=true)
            -> bedtools blacklist filter -> final.bam -> MACS2
```

Duplicates are removed and the blacklist is applied **before** peak calling, not after.

## Blacklist

The default blacklist is downloaded per `--genome`; override it with `--blacklist`:
- **Human (hg38)**: `https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz`
- **Mouse (mm10)**: `https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz`
- Reference: Amemiya et al. 2019 (Scientific Reports, ~1,372 citations)

The hg38 blacklist contains ~900 regions covering ~40 Mb of problematic sequence
including high-signal artifacts, satellite repeats, and assembly gaps. Gzipped BED is
accepted directly.

## Samtools Flag Filtering

For paired-end data the workflow uses `-F 1804`, which removes:
- Bit 4: read unmapped
- Bit 8: mate unmapped
- Bit 256: secondary alignment
- Bit 512: read fails quality checks
- Bit 1024: PCR duplicate (any already flagged by the aligner)

For `--single_end` it uses `-F 1028` (unmapped + duplicate). Neither flag set removes
mitochondrial reads; this workflow does not filter chrM.

## Commands

The workflow runs the equivalent of:

```bash
# Flag and MAPQ filter (paired-end)
samtools view -@ 4 -b -F 1804 -q 30 aligned.bam | \
  samtools sort -@ 4 -o filtered.bam -

# Mark and remove PCR duplicates
picard MarkDuplicates \
  INPUT=filtered.bam \
  OUTPUT=dedup.bam \
  METRICS_FILE=dup_metrics.txt \
  REMOVE_DUPLICATES=true \
  VALIDATION_STRINGENCY=LENIENT

# Remove reads overlapping blacklist regions
bedtools intersect -v -abam dedup.bam -b hg38-blacklist.v2.bed.gz > final.bam
samtools index final.bam

# Record final read count
samtools flagstat final.bam > final_flagstat.txt
```

## Expected Output
- `filtered/<sample>.dup_metrics.txt` -- Picard duplication metrics
- `filtered/<sample>.final.bam` + `.final.bam.bai` -- blacklist-filtered, ready for peak calling
- `filtered/<sample>.final.flagstat.txt` -- final read count after all filtering

The intermediate `filtered.bam` and `dedup.bam` stay in the Nextflow work directory.

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| Duplication rate | <30% (Picard `PERCENT_DUPLICATION`) | Low-input library; consider re-prep |
| Post-filter read count | >=20M TF / >=45M histone | May need deeper sequencing |
| Blacklist overlap | <1% of reads | Expected; higher suggests artifacts |

## Library Complexity Metrics (manual)

NRF, PBC1 and PBC2 are **not computed by this workflow**, and Picard MarkDuplicates does
not report them: its metrics file gives duplicate counts and `PERCENT_DUPLICATION`, while
PBC needs per-position read counts. Compute them from the pre-deduplication BAM if needed:

- **NRF** = distinct read start positions / total reads (Non-Redundant Fraction)
- **PBC1** = positions with exactly 1 read / distinct positions
- **PBC2** = positions with exactly 1 read / positions with exactly 2 reads

```bash
# From the published pre-deduplication BAM (already MAPQ >= 30 from the alignment step)
samtools view -F 1804 -q 30 results/aligned/sample.bam | \
  awk 'BEGIN{OFS="\t"}{print $3, $4, ($2%32>=16 ? "-" : "+")}' | \
  sort | uniq -c | \
  awk '{total += $1; distinct += 1; if ($1 == 1) one += 1; if ($1 == 2) two += 1}
       END {printf "NRF=%.3f PBC1=%.3f PBC2=%.3f\n", distinct/total, one/distinct, one/two}'
```

This is approximate; the ENCODE reference implementation derives the 5' positions with
`bedtools bamtobed` instead of the raw POS field.

These metrics quantify PCR amplification bottleneck severity.
