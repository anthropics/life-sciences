# Filtering, Deduplication, and Spike-in Normalization

CUT&RUN data processing requires standard quality filtering, duplicate removal,
blacklist filtering, AND spike-in normalization for quantitative analysis.

## Quality Filtering

```bash
samtools view -b -h \
    -q 10 \
    -F 1804 \
    -f 2 \
    sample_sorted.bam \
    | samtools sort -@ 4 -o sample_filtered.bam
```

### Filter Parameters

| Flag | Meaning |
|------|---------|
| `-q 10` | MAPQ >= 10 (CUT&RUN uses lower threshold than ChIP-seq) |
| `-F 4` | Remove unmapped |
| `-F 256` | Remove secondary |
| `-F 512` | Remove QC-fail |
| `-F 1024` | Remove duplicates -- a no-op here; duplicates are removed by Picard in the next step |
| `-f 2` | Keep properly paired only |

**Note**: MAPQ 10 instead of 30 for CUT&RUN. The lower threshold retains
more signal because CUT&RUN targets can be in repetitive regions.

## Duplicate Marking

```bash
picard MarkDuplicates \
    INPUT=sample_filtered.bam \
    OUTPUT=sample_dedup.bam \
    METRICS_FILE=sample_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true

samtools index sample_dedup.bam
```

CUT&RUN from low cell numbers may have higher duplication. Accept up to 40%.

## Blacklist Filtering

The workflow filters the BAM against the single file given as `--blacklist`
and stops there: peak files are never filtered, and there is no separate
suspect-list parameter.

```bash
# What the workflow runs, with --blacklist as -b
bedtools intersect \
    -a sample_dedup.bam \
    -b hg38-blacklist.v2.bed \
    -v \
    > sample_final.bam

samtools index sample_final.bam
```

To also exclude the CUT&RUN-specific suspect list (Nordin 2023) -- ~400
regions enriched in CUT&RUN controls that produce false positive peaks,
independent of the ENCODE blacklist -- merge the two files once and pass the
result as `--blacklist`:

```bash
cat hg38-blacklist.v2.bed CUTandRUN.suspectlist.hg38.bed \
    | sort -k1,1 -k2,2n | bedtools merge > combined_blacklist.bed
```

## Spike-in Normalization

### Calculate Scale Factors

The workflow writes one count file per sample
(`spikein/{sample}.spikein_counts.txt`, columns sample and count) and then one
combined `spikein/scale_factors.txt` for the whole run:

```bash
# Per-sample counts, concatenated
cat *.spikein_counts.txt > all_counts.txt

# Scale every sample to the smallest non-zero count (factor = min / count).
# A sample with no spike-in reads is left unscaled (factor 1).
min_count=$(awk -F'\t' '$2 > 0 {print $2}' all_counts.txt | sort -n | head -1)

awk -F'\t' -v min="${min_count:-0}" 'BEGIN {OFS="\t"} {
    factor = ($2 > 0 && min > 0) ? min / $2 : 1
    print $1, $2, factor
}' all_counts.txt > scale_factors.txt
```

### Apply Spike-in Scaling to Signal

The workflow scales the bigWig only, with deepTools:

```bash
# Read scale factor for one sample
scale=$(awk -F'\t' -v s="sample1" '$1==s {print $3}' scale_factors.txt)

bamCoverage \
    --bam sample_final.bam \
    --outFileName sample_normalized.bw \
    --scaleFactor ${scale} \
    --binSize 10 \
    --normalizeUsing None \
    --extendReads \
    --numberOfProcessors 4
```

Without a spike-in index the workflow drops `--scaleFactor` and uses
`--normalizeUsing RPKM` instead.

### Alternative: bedGraph route (not used by the workflow)

```bash
bedtools genomecov \
    -ibam sample_final.bam \
    -bg \
    -pc \
    -scale ${scale} \
    -g hg38.chrom.sizes \
    | sort -k1,1 -k2,2n > sample_normalized.bedGraph

bedGraphToBigWig sample_normalized.bedGraph hg38.chrom.sizes sample_normalized.bw
```

## Generate Fragment BED File

SEACR takes a fragment bedGraph, built from a fragment BED. `bedtools
bamtobed -bedpe` needs mates on adjacent lines, so name-sort the BAM first --
on a coordinate-sorted BAM it warns per read that the mate does not occur next
to it and emits a near-empty BED:

```bash
# Name-sort first
samtools sort -n -@ 2 -o namesorted.bam sample_final.bam

# Keep properly paired fragments on one chromosome and shorter than 1 kb
bedtools bamtobed -bedpe -i namesorted.bam \
    | awk 'BEGIN {OFS="\t"} $1 == $4 && $6 - $2 < 1000 {print $1, $2, $6}' \
    | sort -k1,1 -k2,2n -k3,3n \
    > sample.fragments.bed

# Fragment bedGraph for SEACR -- deliberately unscaled
bedtools genomecov \
    -i sample.fragments.bed \
    -g hg38.chrom.sizes \
    -bg \
    > sample_fragments.bedGraph
```

The workflow does not apply the spike-in factor here: SEACR sees unscaled
fragment coverage, and only the bigWig is calibrated. The BED is published as
`signal/{sample}.fragments.bed`; the bedGraph is an intermediate.

## QC Statistics

```bash
# Final read count
echo "Final reads: $(samtools view -c sample_final.bam)"

# Flagstat (the workflow publishes this as alignment/{sample}.flagstat.txt)
samtools flagstat sample_final.bam > sample_flagstat.txt

# Spike-in fraction; the spike-in BAM is an intermediate, so read the
# published count file instead
genome=$(samtools view -c -F 1804 -f 2 sample_final.bam)
spikein=$(cut -f2 spikein/sample.spikein_counts.txt)
echo "Spike-in fraction: $(echo "scale=4; $spikein / ($genome + $spikein)" | bc)"
```
