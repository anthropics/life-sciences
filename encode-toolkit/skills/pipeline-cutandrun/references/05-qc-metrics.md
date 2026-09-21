# CUT&RUN QC Metrics

Quality assessment for CUT&RUN data includes standard alignment metrics,
CUT&RUN-specific fragment analysis, spike-in validation, and peak quality.

Pass/warn/fail thresholds live in one place: the QC table in `SKILL.md`. The
commands below show how to derive each value, including the ones the workflow
does not compute (spike-in fraction, peak statistics).

## Fragment Size Distribution

The fragment size distribution is the most informative CUT&RUN QC metric:

```bash
# Extract fragment sizes from properly paired reads
samtools view -f 2 -F 1804 sample_final.bam | \
    awk '{if($9 > 0 && $9 < 1000) print $9}' | \
    sort -n | uniq -c | \
    awk '{print $2, $1}' > fragment_sizes.txt
```

The workflow runs exactly this and publishes
`qc/{sample}.fragment_sizes.txt`.

### Expected Patterns by Target

| Target Type | Fragment Pattern | Example |
|-------------|-----------------|---------|
| TF (CTCF, etc.) | Peak <120 bp, some at 150 bp | Sharp sub-nucleosomal |
| Active histone (H3K4me3) | Strong 150 bp peak | Mononucleosomal |
| Repressive histone (H3K27me3) | 150 bp + 300 bp | Mono + dinucleosomal |
| IgG control | Flat distribution | No enrichment pattern |

### Red Flags in Fragment Distribution

- No nucleosomal periodicity: Protocol may have failed
- Only large fragments (>300 bp): Over-digestion or poor tagmentation
- Spike at exact read length: Adapter trimming incomplete
- Identical to IgG: No target enrichment

## Spike-in QC

### Spike-in Fraction

The spike-in BAM is an intermediate; the published per-sample count file is
`spikein/{sample}.spikein_counts.txt` (columns sample and count):

```bash
genome=$(samtools view -c -F 1804 -f 2 results/alignment/sample.filtered.bam)
spikein=$(cut -f2 results/spikein/sample.spikein_counts.txt)
fraction=$(echo "scale=4; $spikein / ($genome + $spikein)" | bc)
echo "Spike-in fraction: $fraction"
```

Interpretation of the fraction (thresholds: see the QC table in `SKILL.md`):
too little spike-in makes the normalization imprecise, too much indicates
poor target enrichment.

### Spike-in Consistency Across Samples

For reliable normalization, spike-in counts should vary across samples
(reflecting different amounts of target material), but not be zero. The
workflow already collects them:

```bash
sort -k2 -n results/spikein/scale_factors.txt
```

Columns: sample, spike-in count, scale factor.

## Alignment Statistics

The workflow publishes `alignment/{sample}.flagstat.txt` and
`alignment/{sample}.dup_metrics.txt`, and the Bowtie2 logs reach
`multiqc/multiqc_report.html`. To recompute:

```bash
samtools flagstat results/alignment/sample.filtered.bam > flagstat.txt
```

Mapping rate, properly paired fraction and duplication rate are judged against
the QC table in `SKILL.md`.

## FRiP (Fraction of Reads in Peaks)

The FRIP process computes this and publishes `qc/{sample}.frip_mqc.tsv`, one
row per peak set called for the sample (SEACR stringent and/or relaxed, and/or
MACS2). Columns: `Peak set` (the peak file), `FRiP`, `reads_in_peaks`,
`total_reads`. MultiQC shows the same table as "Fraction of reads in peaks".

This is what it runs for each peak set:

```bash
total=$(samtools view -c results/alignment/sample.filtered.bam)
in_peaks=$(bedtools intersect \
    -u \
    -a results/alignment/sample.filtered.bam \
    -b results/peaks/sample.seacr.stringent.bed \
    | samtools view -c -)
awk -v a="$in_peaks" -v b="$total" 'BEGIN { printf "FRiP: %.4f\n", a / b }'
```

Both counts are alignments of the filtered BAM with no further flag filter, so
the two mates of a fragment count separately. Counting differently (for
example with `-F 1804 -f 2`) gives a number that does not match the published
one.

CUT&RUN typically has higher FRiP than ChIP-seq because of lower background.
Use the FRiP row of the QC table in `SKILL.md` for the thresholds.

## Peak Count and Size

```bash
# Peak statistics
total_peaks=$(wc -l < results/peaks/sample.seacr.stringent.bed)
echo "Total peaks: $total_peaks"

# Peak size distribution
awk '{print $3-$2}' results/peaks/sample.seacr.stringent.bed | \
    awk '{sum+=$1; n++; a[n]=$1} END {
        asort(a);
        print "Median peak size:", a[int(n/2)];
        print "Mean peak size:", sum/n;
        print "Min:", a[1];
        print "Max:", a[n]
    }'
```

## MultiQC Aggregation

```bash
multiqc \
    --title "CUT&RUN Pipeline QC" \
    --filename multiqc_report \
    --outdir multiqc/ \
    fastqc/ trim_galore/ alignment/ qc/
```

## Summary QC Table

Generate a per-sample summary:

```bash
echo -e "Sample\tReads\tMap_Rate\tDedup_Rate\tSpikein_Frac\tPeaks\tFRiP\tFrag_Peak"
echo -e "${SAMPLE}\t${TOTAL}\t${MAP_RATE}\t${DUP_RATE}\t${SPIKEIN}\t${PEAKS}\t${FRIP}\t${FRAG}"
```
