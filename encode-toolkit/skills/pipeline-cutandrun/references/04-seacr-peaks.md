# SEACR Peak Calling for CUT&RUN

SEACR (Sparse Enrichment Analysis for CUT&RUN) is specifically designed for
the sparse, low-background signal profile of CUT&RUN data. It outperforms
MACS2 on CUT&RUN data because it does not assume a Poisson background model.

## SEACR with IgG Control

When an IgG control is available:

```bash
SEACR_1.3.sh \
    sample_fragments.bedGraph \
    control_IgG_fragments.bedGraph \
    norm \
    stringent \
    sample_seacr
```

### Parameters

| Argument | Value | Description |
|----------|-------|-------------|
| Arg 1 | sample bedGraph | Treatment signal (fragment coverage) |
| Arg 2 | control bedGraph | IgG/no-antibody control signal |
| Arg 3 | `norm` | Normalize control to treatment depth |
| Arg 4 | `stringent` | Stringent mode (conservative peaks) |
| Arg 5 | output prefix | Output file prefix |

## SEACR without Control (Top %)

When no control is available, SEACR uses a numeric threshold:

```bash
SEACR_1.3.sh \
    sample_fragments.bedGraph \
    0.01 \
    non \
    stringent \
    sample_seacr_noctrl
```

`non` is required here: the SEACR v1.3 README states that a numeric threshold
must be paired with `non`. The v1.3 script does not error on `norm` -- it
silently skips normalization -- so the mistake is invisible. The workflow
always passes `non` when there is no `--control`.

The `0.01` value means the top 1% of signal is used as the enrichment threshold.
Adjust based on expected peak count:
- `0.01` (1%): Conservative, fewer peaks
- `0.05` (5%): Moderate
- `0.10` (10%): Permissive, more peaks

## SEACR Stringent vs Relaxed

Run both modes to compare:

```bash
# Stringent: Only peaks that pass both global and local enrichment
SEACR_1.3.sh sample.bedGraph control.bedGraph norm stringent sample_stringent

# Relaxed: Peaks passing global enrichment threshold only
SEACR_1.3.sh sample.bedGraph control.bedGraph norm relaxed sample_relaxed
```

| Mode | Description | Typical Peak Count |
|------|-------------|-------------------|
| Stringent | Global AND local enrichment | 5,000-20,000 |
| Relaxed | Global enrichment only | 10,000-50,000 |

## SEACR Output Format

SEACR v1.3 outputs a 6-column BED-like file:

```
chr1  1000  2000  500.5  8.2  chr1:1450-1480
```

Columns:
1. Chromosome
2. Peak start
3. Peak end
4. Total signal in the peak (AUC)
5. Max signal in the peak
6. Region of maximum signal, as `chr:start-end`

The workflow names these files `{sample}.seacr.stringent.bed` and, with
`--seacr_mode relaxed` or `both`, `{sample}.seacr.relaxed.bed`, published to
`peaks/`.

## Alternative: MACS2 Peak Calling

MACS2 can be used as an alternative or validation:

```bash
macs2 callpeak \
    -t sample_final.bam \
    -c control_IgG.bam \
    -f BAMPE \
    -g hs \
    -n sample_macs2 \
    --nomodel \
    --keep-dup all \
    -q 0.05 \
    --outdir macs2_peaks/
```

### MACS2 Parameters for CUT&RUN

| Parameter | Value | Reason |
|-----------|-------|--------|
| `-f BAMPE` | Paired-end | Use fragment information |
| `--nomodel` | Skip model | CUT&RUN fragments don't follow ChIP model |
| `--keep-dup all` | Keep all | Duplicates already removed |
| `-q 0.05` | FDR 0.05 | Standard threshold |

**Caution**: MACS2 may overcall peaks on CUT&RUN data due to the low
background. Compare with SEACR results and use the intersection for
high-confidence peaks.

## Blacklist + Suspect List Filtering (manual, not run by this workflow)

The workflow applies `--blacklist` to the BAM only and never filters peak
files. To filter the published peaks against both the blacklist and the
CUT&RUN suspect list, run this yourself:

```bash
bedtools intersect \
    -a results/peaks/sample.seacr.stringent.bed \
    -b hg38-blacklist.v2.bed CUTandRUN.suspectlist.hg38.bed \
    -v \
    > sample_peaks_filtered.bed
```

## Peak Overlap Between Callers

When using both SEACR and MACS2, assess concordance:

```bash
# Convert SEACR to 3-column BED
cut -f1-3 results/peaks/sample.seacr.stringent.bed > seacr_peaks.bed
cut -f1-3 results/peaks/sample.macs2_peaks.narrowPeak > macs2_peaks.bed

# Overlap
bedtools intersect -a seacr_peaks.bed -b macs2_peaks.bed -u | wc -l
total_seacr=$(wc -l < seacr_peaks.bed)
total_macs2=$(wc -l < macs2_peaks.bed)
echo "SEACR peaks: $total_seacr"
echo "MACS2 peaks: $total_macs2"
```

Typical overlap: 60-80% of SEACR peaks overlap MACS2 peaks.
High-confidence set: intersection of both callers.

## FRiP Calculation

The FRIP process runs this for every peak set called for a sample and writes
the results to `qc/{sample}.frip_mqc.tsv`:

```bash
total_reads=$(samtools view -c results/alignment/sample.filtered.bam)
reads_in_peaks=$(bedtools intersect \
    -u \
    -a results/alignment/sample.filtered.bam \
    -b results/peaks/sample.seacr.stringent.bed \
    | samtools view -c -)
awk -v a="$reads_in_peaks" -v b="$total_reads" 'BEGIN { printf "FRiP: %.4f\n", a / b }'
```

Both counts are alignments, with no flag filter, so mates count separately.
Judge the result against the QC table in `SKILL.md`; see
`05-qc-metrics.md` for the published file's columns.
