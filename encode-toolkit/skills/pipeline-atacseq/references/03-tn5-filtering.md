# Stage 3: Filtering, Tn5 Shift, and Fragment Selection

## Tools
- **Samtools 1.19** (image version): Mitochondrial read removal, sorting, indexing
- **Picard MarkDuplicates 3.1.1** (image version, Java 17): PCR duplicate removal
- **deeptools `alignmentSieve` 3.5.5** (image version): Tn5 shift and fragment size selection
- **bedtools 2.31.0** (image version): Blacklist filtering

## Order of operations in the workflow

```
aligned.bam -> mito removal -> Picard MarkDuplicates (REMOVE_DUPLICATES=true)
            -> alignmentSieve --ATACshift -> bedtools blacklist filter
            -> final.bam -> alignmentSieve size selection (NFR, mono-nucleosome)
```

The Tn5 shift runs **after** duplicate removal, not before. Fragment size selection is the
last step, and peaks are called on the NFR BAM.

## Tn5 Transposase Offset Correction

The Tn5 transposase creates a 9-bp target site duplication during insertion. To
accurately represent the cut site, reads must be shifted:
- **Forward strand (+)**: shift +4 bp
- **Reverse strand (-)**: shift -5 bp

This correction is critical for motif footprinting and accurate cut-site analysis.

## Commands

The workflow runs the equivalent of:

```bash
# Step 1: Record the per-contig read counts, then remove chrM reads
samtools idxstats aligned.bam > sample.idxstats.txt
samtools view -@ 4 -b aligned.bam $(samtools idxstats aligned.bam | \
  awk '$1 != "chrM" && $1 != "*" {print $1}' | tr '\n' ' ') > no_mito.bam

# Step 2: Mark and remove PCR duplicates
samtools sort -@ 4 -o sorted.bam no_mito.bam
picard MarkDuplicates \
  INPUT=sorted.bam OUTPUT=dedup.bam \
  METRICS_FILE=dup_metrics.txt \
  REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
samtools index dedup.bam

# Step 3: Apply Tn5 shift (+4/-5)
alignmentSieve --bam dedup.bam --outFile shifted_unsorted.bam \
  --ATACshift --numberOfProcessors 4
samtools sort -@ 4 -o shifted.bam shifted_unsorted.bam
samtools index shifted.bam

# Step 4: Remove reads overlapping blacklist regions
bedtools intersect -v -abam shifted.bam -b hg38-blacklist.v2.bed.gz > final.bam
samtools index final.bam
samtools flagstat final.bam > final.flagstat.txt

# Step 5: Separate nucleosome-free and mono-nucleosomal fragments
alignmentSieve --bam final.bam --outFile nfr_unsorted.bam \
  --maxFragmentLength 150 --numberOfProcessors 4
samtools sort -@ 4 -o nfr.bam nfr_unsorted.bam && samtools index nfr.bam

alignmentSieve --bam final.bam --outFile mononuc_unsorted.bam \
  --minFragmentLength 150 --maxFragmentLength 300 --numberOfProcessors 4
samtools sort -@ 4 -o mononuc.bam mononuc_unsorted.bam
```

The fragment length boundary (150) is `--nfr_max`, and the mitochondrial contig name is
`--mito_name`. There is no `samtools view -F 1804` step in this workflow: flag filtering
(`-q 30 -f 2`) happens once, during alignment.

## Expected Output
- `qc/<sample>.idxstats.txt` -- `samtools idxstats` of the BAM **before** mitochondrial
  reads are removed: one row per contig with name, length, mapped and unmapped counts. The
  mitochondrial fraction is the mapped count on the `--mito_name` row divided by the sum of
  the mapped column; MultiQC's samtools module reads the same file and reports it. The step
  also works when the genome has no such contig -- nothing is removed. No `bc` is involved.
- `filtered/<sample>.dup_metrics.txt` -- Picard duplication metrics
- `filtered/shifted/<sample>.shifted.bam` + `.bai` -- Tn5-corrected, pre-blacklist
- `filtered/<sample>.final.bam` + `.final.bam.bai` -- blacklist-filtered (all fragments)
- `filtered/<sample>.final.flagstat.txt` -- read count after all filtering
- `filtered/nfr/<sample>.nfr.bam` + `.bai` -- nucleosome-free fragments
- `filtered/nfr/<sample>.mononuc.bam` -- mono-nucleosome fragments (no index is written)

Both size-selected BAMs land in `filtered/nfr/`.

## Fragment Size Classes

| Class | Size (bp) | Use | Produced here? |
|-------|-----------|-----|----------------|
| NFR (nucleosome-free) | <150 (`--nfr_max`) | Peak calling, TF footprinting | yes |
| Mono-nucleosome | 150-300 | Nucleosome positioning | yes |
| Di-nucleosome | 300-500 | Chromatin architecture | no |

## QC Checkpoints

| Check | Threshold | Action if Failed |
|-------|-----------|------------------|
| Mitochondrial fraction (`qc/<sample>.idxstats.txt`) | <20% (ideal <5%) | Optimize cell lysis |
| Duplication rate | <30% (Picard `PERCENT_DUPLICATION`) | Low complexity library |
| NFR fraction | >40% of fragments <150bp (manual) | Check transposition efficiency |
| Post-filter reads | >=25M | May need deeper sequencing |

## Notes
- `alignmentSieve --ATACshift` applies the +4/-5 offset automatically. It emits an unsorted
  BAM, so every `alignmentSieve` call is followed by `samtools sort`.
- NRF and PBC are not computed by this workflow; see `references/05-qc-metrics.md`.
- Blacklist: Amemiya et al. 2019 (hg38-blacklist.v2.bed.gz, ~900 regions, ~40 Mb). It is
  applied to the BAM, so peaks are already blacklist-clean.
