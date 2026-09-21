# TF Footprinting Analysis for DNase-seq

Transcription factor footprinting detects protein-DNA interactions from
DNase-seq cleavage patterns. Bound TFs protect DNA from cleavage, creating
a "footprint" -- a local depression in the DNase-seq signal.

## Prerequisites

Reliable footprinting requires:
- **Deep sequencing**: >100 million uniquely mapped reads
- **Paired-end data**: Better resolution than single-end
- **High SPOT score**: >0.4 (clean signal)
- **Peak calls**: DHS regions from Hotspot2
- **An RGT data directory**: HINT reads genome sequence and annotation from
  `$RGTDATA`; the workflow stages the directory passed with `--rgt_data`

## HINT Footprinting (DNase-seq mode)

HINT (from the Regulatory Genomics Toolbox) works for both DNase-seq and
ATAC-seq data. This is what the workflow runs:

```bash
rgt-hint footprinting \
    --dnase-seq \
    --paired-end \
    --organism hg38 \
    --output-location footprints/ \
    --output-prefix sample \
    sample.filtered.bam \
    sample.peaks.narrowPeak
```

It writes `footprints/sample.bed`, which the workflow publishes as
`results/footprints/sample.footprints.bed`. No other footprinting output is
produced.

### Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--dnase-seq` | Flag | Use DNase-seq cleavage model (not ATAC-seq) |
| `--paired-end` | Flag | Use paired-end fragment information |
| `--organism` | hg38 | Genome build for bias correction; must be set up in the RGT data directory |
| `--output-prefix` | sample | Base name of the output BED |

**Important**: Use `--dnase-seq` for DNase-seq data and `--atac-seq` for
ATAC-seq data. They have different cleavage bias models, and using the wrong
one silently produces wrong footprints.

## Motif Matching in Footprints (manual, optional)

Motif matching is not part of the workflow. To run it afterwards, note that
`rgt-motifanalysis matching --motif-dbs` takes **directories of `.pwm` files**
in the pre-2016 JASPAR format, not a `.meme` file. The RGT data directory ships
such directories under `motifs/`:

```bash
export RGTDATA=/ref/rgtdata

rgt-motifanalysis matching \
    --organism hg38 \
    --input-files results/footprints/sample.footprints.bed \
    --output-location motif_matches/ \
    --motif-dbs /ref/rgtdata/motifs/jaspar_vertebrates
```

To use a motif collection RGT does not ship, convert it to `.pwm` files in
their own directory first with RGT's `createPwm.py`, then point `--motif-dbs`
at that directory.

## Wellington Footprinting (alternative, manual)

Wellington (pyDNase) uses a different statistical approach. It is **not** in the
pipeline image or in `dnaseseq-env.yml`; install it separately with
`pip install pyDNase`.

```bash
wellington_footprints.py \
    -A \
    -p 20 \
    -fdrlimit 0.01 \
    sample.peaks.narrowPeak \
    sample.filtered.bam \
    wellington_out/
```

### HINT vs Wellington Comparison

| Feature | HINT | Wellington |
|---------|------|-----------|
| Bias correction | Sequence-specific | Position-based |
| Speed | Moderate | Fast |
| Sensitivity | Higher | More conservative |
| DNase + ATAC | Both | Both |
| Active development | Yes | Limited |
| In the pipeline image | Yes (RGT 1.0.2) | No (`pip install pyDNase`) |

## Footprint Quality Assessment

None of the checks below are computed by the workflow; run them on the
published footprint BED.

### Per-Motif Footprint Depth

```bash
# Calculate average footprint score at known CTCF sites
bedtools intersect \
    -a CTCF_motif_sites.bed \
    -b results/footprints/sample.footprints.bed \
    -wa -wb \
    | awk '{print $NF}' \
    | awk '{sum+=$1; n++} END {print "Mean CTCF footprint score:", sum/n}'
```

### Footprint vs Background Signal Ratio

Good footprints show:
- Clear signal depression at the motif center
- Flanking shoulders of higher cleavage
- Depth-to-flank ratio > 1.5

### Expected Footprint Counts

| Sequencing Depth | Expected Footprints |
|------------------|---------------------|
| 50M reads | Unreliable |
| 100M reads | 50,000-100,000 |
| 200M reads | 100,000-200,000 |
| 500M reads | 200,000-400,000 |

## Aggregate Footprint Visualization (manual)

Generate aggregate footprint profiles across all instances of a motif. This
needs the motif matches from the step above:

```bash
rgt-hint differential \
    --organism hg38 \
    --bc \
    --nc 8 \
    --mpbs-files motif_matches/sample_mpbs.bed \
    --reads-files sample.filtered.bam \
    --conditions sample \
    --output-location diff_footprints/
```

This produces per-motif aggregate profiles showing the average cleavage
pattern across all binding sites, which is more robust than individual
footprint calls.

## Comparison with Vierstra 2020 Reference Map

The Vierstra et al. 2020 reference map provides a global catalog of human
TF footprints from 243 DNase-seq datasets. Use it to:

1. Validate your footprint calls
2. Compare tissue-specific footprinting
3. Identify novel TF binding events

```bash
# Download Vierstra reference footprints
# Available at: https://www.vierstra.org/resources/dgf

# Compare overlap
bedtools intersect \
    -a results/footprints/sample.footprints.bed \
    -b vierstra_consensus_footprints.bed \
    -u | wc -l
```
