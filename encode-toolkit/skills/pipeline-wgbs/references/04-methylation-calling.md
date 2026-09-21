# Methylation Calling with MethylDackel

MethylDackel (formerly PileOMeth) extracts per-CpG methylation levels from
bisulfite-aligned BAM files. This workflow then converts its bedGraph output to the
ENCODE bedMethyl layout (https://www.encodeproject.org/data-standards/wgbs/).

## M-bias Assessment (Run First)

The workflow runs `MethylDackel mbias` on every sample before extraction:

```bash
MethylDackel mbias \
    --CHG --CHH \
    /ref/genome/genome.fa \
    sample.sorted.bam \
    sample_mbias > sample_mbias_report.txt 2>&1
```

This produces SVG plots showing methylation level by read position, and the redirected
text output carries MethylDackel's suggested inclusion bounds. Both land in
`bismark/mbias/`. Look for:
- Elevated methylation at the 5' end of read 2 (end-repair artifact). The workflow already
  removes this at trimming with `--clip_R2 10`, so it should be gone
- Irregular methylation at read ends (adapter contamination)

## Methylation Extraction

This is exactly what the workflow runs (`--mergeContext` is dropped when
`--merge_context false`):

```bash
MethylDackel extract \
    --mergeContext \
    --CHG --CHH \
    --opref sample \
    /ref/genome/genome.fa \
    sample.sorted.bam
```

### Key Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--mergeContext` | Enabled by `--merge_context true` | Merge forward/reverse strand CpG data into one record |
| `--CHG --CHH` | Enabled | Also extract non-CpG methylation (useful for ESCs, neurons) |
| `--opref sample` | Output prefix | Names the three bedGraphs |

### Parameters the Workflow Deliberately Does Not Pass

| Parameter | Why not |
|-----------|---------|
| `--minDepth` | Extraction is unfiltered so the coverage statistics can see every covered site. The `--min_coverage` cut is applied when the bedMethyl files are written |
| `--maxDepth` | Not an option in MethylDackel 0.6.1. `getopt_long` rejects the long form and `extract` exits with a usage error. Do not pass it |
| `--nOT` / `--nOB` | The 5' end-repair bias of read 2 is removed at trimming (`--clip_R2 10`), so no read positions are excluded here |

### How `--OT`/`--nOT` Arguments Are Ordered

Getting this wrong silently trims the wrong end. From the MethylDackel 0.6.1 help,
`--OT A,B,C,D` means "include calls at positions from A through B on read #1 and C
through D on read #2". The `--nOT a,b,c,d` form is the complement and *excludes*:

| Slot | Meaning |
|------|---------|
| 1 | bases from the start (5') of read 1 |
| 2 | bases from the end (3') of read 1 |
| 3 | bases from the start (5') of read 2 |
| 4 | bases from the end (3') of read 2 |

So `--nOT 0,0,0,10` trims read 2's **3'** end, not its 5' end; read 2's 5' end is
`--nOT 0,0,10,0`. `--nOB 0,10,0,0` trims read 1's 3' end. If the M-bias plots still show
a problem, re-extract by hand with the bounds the mbias report suggests.

### Context-Specific Output Files

MethylDackel produces one bedGraph per context, named from `--opref` with an underscore.
They are published to `bismark/methylation/` and hold every covered site:

- `<sample>_CpG.bedGraph` -- CpG methylation (primary)
- `<sample>_CHG.bedGraph` -- CHG methylation (non-CpG)
- `<sample>_CHH.bedGraph` -- CHH methylation (non-CpG)

Columns are chrom, start, end, methylation percentage, methylated read count,
unmethylated read count, and the file starts with a `track` header line.

## Convert to bedMethyl Format

The workflow converts each context with this awk, then sorts, compresses and indexes:

```bash
for context in CpG CHG CHH; do
    awk -v min_cov=5 'BEGIN {OFS="\t"} !/^track/ {
        cov = $5 + $6
        if (cov == 0 || cov < min_cov) next
        score = (cov > 1000) ? 1000 : cov
        print $1, $2, $3, ".", score, ".", $2, $3, "0,0,0", cov, int(($5 / cov) * 100 + 0.5)
    }' sample_${context}.bedGraph \
        | sort -k1,1 -k2,2n \
        | bgzip > sample.${context}.bedMethyl.gz

    tabix -p bed sample.${context}.bedMethyl.gz
done
```

Four details matter and are easy to get wrong:

- **`!/^track/`** skips MethylDackel's header line. Without it the header is parsed as data.
- **Column 5 (score)** is the read count capped at 1000, per the ENCODE definition — not
  the methylation percentage scaled by 10.
- **Column 6 (strand)** is `.`. The MethylDackel bedGraph carries no strand, and with
  `--mergeContext` a record covers both strands of the CpG anyway.
- **Column 11** is an integer percentage (rounded), and column 10 is the coverage.

Sites below `min_cov` are dropped, which is the only place `--min_coverage` is applied.

## Per-Chromosome Extraction (Parallel, not run by this workflow)

For very large BAMs you can parallelize extraction by chromosome by hand:

```bash
for chr in $(samtools idxstats sample.sorted.bam | cut -f1 | grep -v '*'); do
    MethylDackel extract \
        --mergeContext \
        --opref "perchr/${chr}" \
        -r "${chr}" \
        /ref/genome/genome.fa \
        sample.sorted.bam &
done
wait

# Concatenate results, dropping the per-file track headers
cat perchr/*_CpG.bedGraph | awk '!/^track/' | sort -k1,1 -k2,2n > sample_CpG.bedGraph
```

## Global Methylation Summary (manual)

Compute genome-wide methylation statistics from a published bedGraph:

```bash
awk '!/^track/ {
    meth += $5; unmeth += $6; n++
} END {
    total = meth + unmeth;
    print "Covered sites:", n;
    print "Mean methylation:", (meth/total)*100 "%";
    print "Methylated reads:", meth;
    print "Unmethylated reads:", unmeth
}' sample_CpG.bedGraph
```

Expected values for mammalian somatic tissue:
- Global CpG methylation: 70-85%
- CpG islands: 5-15% (mostly unmethylated)
- Gene bodies: 60-80%
- Intergenic: 75-90%
