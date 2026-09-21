# Stage 4: Signal Track Generation

## Tools
- **STAR**: Generates strand-specific bedGraph output during alignment
- **bedGraphToBigWig**: Converts bedGraph to bigWig (UCSC Kent tools)

## Overview

RNA-seq signal tracks display read coverage across the genome for visualization in
genome browsers (UCSC, IGV, WashU Epigenome Browser). For stranded RNA-seq, separate
plus-strand and minus-strand tracks are generated, enabling gene-level visualization
of sense and antisense transcription.

## STAR bedGraph Output

STAR generates bedGraph files directly during alignment. The workflow passes
`--outWigType bedGraph` always, and `--outWigStrand` follows `--strandedness`:

| `--strandedness` | STAR receives | bedGraph files written |
|------------------|---------------|------------------------|
| `reverse`, `forward` | `--outWigStrand Stranded` | `Signal.{Unique,UniqueMultiple}.str1.out.bg` and `...str2.out.bg` |
| `none` | `--outWigStrand Unstranded` | `Signal.{Unique,UniqueMultiple}.str1.out.bg` only |

## What str1 and str2 Mean

STAR assigns a fragment to str1 when **read 1 maps to the + strand** of the genome
(`signalFromBAM.cpp`: `iStrand = ((flag&0x10)>0) == ((flag&0x80)==0)`). str1 is therefore
not "the plus strand" in the transcript sense — which transcript strand it represents
depends on the library:

| `--strandedness` | Library | str1 holds | str2 holds |
|------------------|---------|-----------|-----------|
| `reverse` (dUTP, ENCODE standard) | read 1 is antisense to the transcript | minus-strand transcription | plus-strand transcription |
| `forward` (directional ligation) | read 1 is sense to the transcript | plus-strand transcription | minus-strand transcription |
| `none` | unstranded | all fragments, one track | — |

The workflow applies exactly this mapping, matching ENCODE's `STAR_RSEM.sh`
(`str[1]="-"; str[2]="+"` for the dUTP default), and names the outputs accordingly:

| `--strandedness` | Published tracks |
|------------------|------------------|
| `reverse` | str1 -> `signal/<sample>_minus.bw`, str2 -> `signal/<sample>_plus.bw` |
| `forward` | str1 -> `signal/<sample>_plus.bw`, str2 -> `signal/<sample>_minus.bw` |
| `none` | str1 -> `signal/<sample>_unstranded.bw` |

## Unique vs UniqueMultiple

`Signal.UniqueMultiple.*` includes multi-mapped reads (ENCODE standard) and is what the
workflow converts to bigWig. `Signal.Unique.*` counts only uniquely mapped reads and is a
more conservative estimate; it is published as
`star/<sample>.Signal.Unique.str*.out.bg` but is not converted, so run the commands below
on it yourself if you want that track as a bigWig.

## bedGraph to bigWig Conversion

This is what the workflow runs for a `reverse`-stranded sample:

```bash
# Sort bedGraph (required by bedGraphToBigWig)
sort -k1,1 -k2,2n sample.Signal.UniqueMultiple.str1.out.bg > minus_sorted.bg
sort -k1,1 -k2,2n sample.Signal.UniqueMultiple.str2.out.bg > plus_sorted.bg

# Convert to bigWig
bedGraphToBigWig minus_sorted.bg chrom.sizes sample_minus.bw
bedGraphToBigWig plus_sorted.bg chrom.sizes sample_plus.bw
```

Swap the two strand labels for a `forward`-stranded library.

## Chromosome Sizes File

The workflow uses `--chrom_sizes`, which defaults to `chrNameLength.txt` inside the STAR
index directory. To build one by hand:

```bash
# From the aligned BAM header
samtools view -H sample.Aligned.sortedByCoord.out.bam | \
  grep '@SQ' | awk '{print $2"\t"$3}' | \
  sed 's/SN://;s/LN://' > chrom.sizes

# Or fetch from UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```

## Signal Normalization

STAR bedGraph output is raw read counts per position, and the workflow publishes the
bigWigs unnormalized. For cross-sample comparison, normalize them yourself afterwards:

| Method | Description | When to Use |
|--------|-------------|-------------|
| Raw | Unnormalized read counts (what this workflow writes) | Single-sample visualization |
| RPM | Reads per million mapped | Cross-sample comparison (simple) |
| RPKM | RPM per kilobase | Length-normalized (rarely needed for signal) |

ENCODE distributes both raw and RPM-normalized bigWig files. For publication, RPM is
standard for cross-sample comparison in browser screenshots.

## Manual Alternative: bamCoverage (deepTools, not run by this workflow)

deepTools is in neither the container image nor `rnaseq-env.yml`, so this needs a
separate install. It gives more control over normalization and strand selection:

```bash
bamCoverage -b sample.Aligned.sortedByCoord.out.bam \
  -o sample.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --filterRNAstrand forward \
  --numberOfProcessors 8
```

## Expected Output
- `signal/<sample>_plus.bw` -- plus-strand transcription (stranded runs)
- `signal/<sample>_minus.bw` -- minus-strand transcription (stranded runs)
- `signal/<sample>_unstranded.bw` -- single track when `--strandedness none`

## Notes
- bedGraphToBigWig requires sorted input and a chromosome sizes file.
- For stranded data, always keep separate plus/minus tracks. Combining them
  loses strand information and confounds sense/antisense transcription.
- A swapped plus/minus pair is easy to miss: check a gene you know the orientation of
  in a browser before publishing tracks.
- bigWig files are typically 50-200 MB each, much smaller than BAM files.
