# WGBS QC Metrics and Conversion Rate Assessment

Quality control for WGBS requires bisulfite-specific metrics beyond
standard alignment QC. The most critical metric, the bisulfite conversion rate, is not
computed by this workflow — the sections below say which commands you have to run
yourself.

## What the Workflow Produces

| Output | Contents |
|--------|----------|
| `bismark/alignments/*_PE_report.txt` | Mapping efficiency, and the percentage of C methylated in CpG/CHG/CHH context |
| `bismark/dedup_reports/*.deduplication_report.txt` | Duplicate count and rate |
| `bismark/mbias/<sample>_mbias_*.svg`, `<sample>_mbias_report.txt` | M-bias plots and MethylDackel's suggested inclusion bounds |
| `coverage/<sample>.coverage_stats.txt` | Covered CpGs, their mean coverage, and the fraction reaching >=5x, >=10x and `--min_coverage`x |
| `multiqc/multiqc_report.html` | Aggregated report |

## Coverage Statistics

`COVERAGE_STATS` parses the unfiltered `<sample>_CpG.bedGraph`, so it sees every covered
site rather than only the ones that pass `--min_coverage`:

```bash
awk -v min_cov=5 '!/^track/ {
    cov = $5 + $6; sum += cov; n++;
    if (cov >= 5)       c5++;
    if (cov >= 10)      c10++;
    if (cov >= min_cov) cmin++
} END {
    if (n == 0) { print "Covered CpGs: 0"; exit }
    printf "Covered CpGs (>=1x): %d\n", n;
    printf "Mean coverage of covered CpGs: %.1f\n", sum/n;
    printf "Covered CpGs >=5x: %d (%.1f%%)\n", c5, c5/n*100;
    printf "Covered CpGs >=10x: %d (%.1f%%)\n", c10, c10/n*100;
    printf "Covered CpGs >=%dx (--min_coverage, kept in bedMethyl): %d (%.1f%%)\n", min_cov, cmin, cmin/n*100
}' sample_CpG.bedGraph > sample.coverage_stats.txt
```

`min_cov` is `--min_coverage`, so the five-line file looks like this at the default of 5:

```
Covered CpGs (>=1x): 27184023
Mean coverage of covered CpGs: 12.4
Covered CpGs >=5x: 22903511 (84.3%)
Covered CpGs >=10x: 16992841 (62.5%)
Covered CpGs >=5x (--min_coverage, kept in bedMethyl): 22903511 (84.3%)
```

Three things to keep straight when reporting these numbers:

- The 5x and 10x thresholds are fixed in the workflow. They do not follow `--min_coverage`.
- The last line is the one that does, and it is the count that reached the bedMethyl
  files. At the default it duplicates the >=5x line; at any other `--min_coverage` it does
  not.
- The percentages are of *covered* CpGs. A CpG that received zero reads is not in the
  bedGraph and is not counted, so this is not the genome-wide CpG completeness.

With `--merge_context false` the records are per cytosine, not per CpG, and the counts
roughly double.

## Bisulfite Conversion Rate (Lambda Spike-in, not run by this workflow)

Lambda phage DNA is fully unmethylated. Any methylation detected on lambda
represents incomplete bisulfite conversion. Run this yourself, against the trimmed reads,
if a spike-in was included:

```bash
# Align to lambda genome
bismark \
    --genome /ref/lambda/ \
    --bowtie2 \
    --parallel 2 \
    -1 sample_R1_val_1.fq.gz \
    -2 sample_R2_val_2.fq.gz \
    --output_dir lambda_out/ \
    --unmapped

# Extract methylation from lambda alignments
MethylDackel extract \
    --mergeContext \
    --minDepth 1 \
    /ref/lambda/genome.fa \
    lambda_out/sample_pe.bam

# Calculate conversion rate
awk '!/^track/ {meth+=$5; unmeth+=$6} END {
    total=meth+unmeth;
    conv=(unmeth/total)*100;
    print "Conversion rate: " conv "%";
    print "Unconverted (false methylation): " (meth/total)*100 "%"
}' lambda_CpG.bedGraph
```

### Conversion Rate Thresholds

| Rate | Status | Action |
|------|--------|--------|
| >99.5% | Excellent | Proceed |
| 99.0-99.5% | Acceptable | Proceed with note |
| 98.0-99.0% | Warning | May inflate methylation estimates |
| <98.0% | Fail | Do NOT use this library |

## Non-CpG Methylation as Conversion Proxy

If no spike-in is available, use CHH methylation as a proxy. The Bismark alignment report
already carries "C methylated in CHH context"; the same number can be recomputed from the
published CHH bedGraph:

```bash
awk '!/^track/ {meth+=$5; unmeth+=$6} END {
    print "CHH methylation: " (meth/(meth+unmeth))*100 "%"
}' sample_CHH.bedGraph
```

In somatic tissue, CHH methylation should be <1%. Higher values suggest
incomplete conversion. Exception: embryonic stem cells and neurons can have
genuine non-CpG methylation (2-5%), which is why this is a proxy rather than a
measurement.

## Additional Manual Coverage Checks

```bash
# Genome-wide coverage distribution from the published BAM
samtools depth -a bismark/alignments/sample.sorted.bam | \
    awk '{cov[$3]++} END {for (c in cov) print c, cov[c]}' | \
    sort -k1,1n > coverage_distribution.txt

# Mean and median coverage (needs gawk for asort)
samtools depth -a bismark/alignments/sample.sorted.bam | \
    awk '{sum+=$3; n++; a[n]=$3; if ($3>=5) sum5++} END {
        asort(a);
        print "Mean:", sum/n;
        print "Median:", a[int(n/2)];
        print "Total bases:", n;
        print "Bases >=5x:", sum5/n*100 "%"
    }'
```

## Mapping Statistics

The mapping efficiency in the Bismark `*_PE_report.txt` is the primary number. For a
flag-level breakdown of the final BAM:

```bash
samtools flagstat bismark/alignments/sample.sorted.bam > flagstat.txt
```

Key values to extract:
- Total reads in the final (deduplicated) BAM
- Mapped reads (Bismark mapping efficiency >70% expected)
- Properly paired (expect >95% of mapped)
- Duplication rate (from the deduplication report)

## MultiQC Report

The workflow runs:

```bash
multiqc --title "ENCODE WGBS Pipeline" --filename multiqc_report --force .
```

`--filename` is required: without it MultiQC 1.21 derives the file name from `--title`
and writes `ENCODE-WGBS-Pipeline_multiqc_report.html`, which would not match the declared
output.

The inputs collected are the raw-read FastQC reports, the Trim Galore trimming reports,
the FastQC reports for the trimmed reads, the Bismark `*_PE_report.txt` files and the
`*.deduplication_report.txt` files, so the duplication rate appears in the report as well
as in `bismark/dedup_reports/`. With `--skip_dedup` there are no deduplication reports to
collect. Picard metrics do not exist for this pipeline.

## Summary QC Table Format

Generate a per-sample summary for reporting:

```bash
echo -e "Sample\tTotal_Reads\tMapping_Rate\tDedup_Rate\tConversion\tMean_CpG_Cov\tCpGs_5x"
echo -e "${SAMPLE}\t${TOTAL}\t${MAP_RATE}\t${DEDUP_RATE}\t${CONV_RATE}\t${MEAN_COV}\t${CPGS_5X}"
```

`CONV_RATE` has to come from the manual check above; everything else is in the published
reports.
