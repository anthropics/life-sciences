# Hotspot2 DHS Calling and Signal Generation

Hotspot2 is the ENCODE-standard peak caller for DNase-seq data. It identifies
DNase I Hypersensitive Sites (DHSs) using a local tag density model that
accounts for mappability variation across the genome.

## Hotspot2 Execution

`hotspot2.sh` takes options followed by two positional arguments: the input BAM
and the output directory. There is no `-s` (input) or `-o` (output) option.
`-c` expects the chromosome sizes as a sorted BED file (chromosome, 0, length),
not the two-column `.chrom.sizes` file, so convert it first. This is exactly what
`scripts/main.nf` runs:

```bash
# hotspot2.sh wants chromosome sizes as a sorted BED file with column 2 set to 0
awk 'BEGIN {OFS="\t"} {print $1, 0, $2}' hg38.chrom.sizes | sort-bed - > chrom_sizes.bed

hotspot2.sh \
    -c chrom_sizes.bed \
    -C /ref/hotspot2/hg38.center_sites.n100.starch \
    -M /ref/hotspot2/hg38.mappable_only.bed \
    -f 0.05 \
    -F 0.05 \
    sample.filtered.bam \
    hotspot2_out
```

### Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-c` | `chrom_sizes.bed` | Required. Chromosome sizes as a sorted BED file: chromosome, 0, length |
| `-C` | `center_sites.starch` | Required. Center sites, made once per genome by `extractCenterSites.sh` |
| `-M` | `mappable_regions.bed` | Optional. The mappable-regions BED the center sites were derived from (read-length specific) |
| `-f` | 0.05 | Hotspot FDR threshold. Every output file is named after this value |
| `-F` | 0.05 | Site-call FDR threshold. Must be >= `-f`, or `hotspot2.sh` aborts |
| `-p` | default | Peak definition (peak shape), not a protocol name. The workflow leaves it at the default |
| positional 1 | `in.bam` | Input BAM |
| positional 2 | `outdir` | Output directory |

`-F` may not be stricter than `-f`. Because `--fdr` is user-supplied, `main.nf`
passes `max(--fdr, 0.05)` for `-F` so that a loose `--fdr` cannot abort the run.
If you call `hotspot2.sh` by hand with a variable FDR, do the same:

```bash
SITECALL_FDR=$(echo "$FDR 0.05" | awk '{print ($1 > $2) ? $1 : $2}')
hotspot2.sh -c chrom_sizes.bed -C center_sites.starch -f "$FDR" -F "$SITECALL_FDR" in.bam outdir
```

### Center Sites and Mappability

Hotspot2 needs a center-sites archive, which is derived from a mappable-regions
BED. Both are genome-build and read-length specific. There is no
`hotspot2-mappability` program; the script that ships with Hotspot2 is
`extractCenterSites.sh`, and it is on the PATH inside the pipeline image:

```bash
extractCenterSites.sh \
    -c chrom_sizes.bed \
    -M hg38.mappable_only.bed \
    -o hg38.center_sites.n100.starch
```

Mappable-regions BEDs are produced separately (for example with Umap/Bismap or
the ENCODE-provided files) for the read length of the library. Pass the same
file to `hotspot2.sh -M` that was used to build the center sites.

## Output Files

`hotspot2.sh` names every output after the input BAM's basename and the `-f`
threshold, never `-F`. For `sample.filtered.bam` with `-f 0.05`:

| File | Format | Description |
|------|--------|-------------|
| `sample.filtered.hotspots.fdr0.05.starch` | starch | Hotspot regions passing the FDR cutoff |
| `sample.filtered.peaks.fdr0.05.starch` | starch | Individual peaks within hotspots |
| `sample.filtered.peaks.fdr0.05.narrowpeaks.starch` | starch | The same peaks in narrowPeak columns |
| `sample.filtered.SPOT.fdr0.05.txt` | text | SPOT score |
| `sample.filtered.allcalls.starch` | starch | All site calls before FDR filtering |
| `sample.filtered.density.starch` | starch | Per-base cleavage density |
| `sample.filtered.density.bw` | bigWig | The same cleavage density as a bigWig |

Every `.starch` file is a compressed BEDOPS archive, not text. Run `unstarch`
before feeding it to `bedtools`, `awk`, or anything else that expects BED.

## Convert to BED and narrowPeak

```bash
unstarch hotspot2_out/sample.filtered.hotspots.fdr0.05.starch > sample.hotspots.fdr0.05.bed
unstarch hotspot2_out/sample.filtered.peaks.fdr0.05.narrowpeaks.starch > sample.peaks.narrowPeak
unstarch hotspot2_out/sample.filtered.allcalls.starch > sample.allcalls.bed
cp hotspot2_out/sample.filtered.SPOT.fdr0.05.txt sample.SPOT.txt
```

These four files are what the workflow publishes to `results/hotspots/`.

The workflow filters the BAM against the blacklist before Hotspot2 runs, so the
peaks are already blacklist-free. A peak-level intersect is only needed if you
call hotspots on a BAM that was not filtered, or if you want to apply an
additional list:

```bash
bedtools intersect -a sample.peaks.narrowPeak -b /ref/hg38-blacklist.v2.bed -v \
    > sample.DHS.narrowPeak
```

## SPOT Score

The SPOT score (Signal Portion of Tags) is computed by Hotspot2:

```bash
cat sample.SPOT.txt
```

If computing manually from the converted BED:
```bash
total_tags=$(samtools view -c sample.filtered.bam)
tags_in_hotspots=$(bedtools intersect \
    -a sample.filtered.bam \
    -b sample.hotspots.fdr0.05.bed \
    -u -bed | wc -l)
echo "SPOT score: $(echo "scale=4; $tags_in_hotspots / $total_tags" | bc)"
```

## Signal Track Generation

Hotspot2 writes its own cleavage-density bigWig (`*.density.bw`, above). The
workflow does not publish it; it builds a separate read-per-million fragment
coverage track with bedtools and publishes that as
`results/hotspots/<sample>.density.bw`:

```bash
# RPM scale factor from the filtered, properly paired reads
total=$(samtools view -c -F 1804 -f 2 sample.filtered.bam)
scale=$(echo "scale=10; 1000000 / $total" | bc)

bedtools genomecov \
    -ibam sample.filtered.bam \
    -bg \
    -pc \
    -g /ref/hg38.chrom.sizes \
    | awk -v s=$scale 'BEGIN{OFS="\t"} {$4=$4*s; print}' \
    | sort -k1,1 -k2,2n \
    > sample_rpm.bedGraph

bedGraphToBigWig sample_rpm.bedGraph /ref/hg38.chrom.sizes sample.density.bw
```

Use the Hotspot2 density track instead when you want per-base cut counts rather
than fragment coverage.

## Peak Annotation

Annotate DHSs with genomic features (a manual downstream step; the workflow
stops at peak calling):

```bash
# Count peaks by category
total=$(wc -l < sample.peaks.narrowPeak)
echo "Total DHSs: $total"

# Overlap with gene promoters (TSS +/- 2kb)
promoter=$(bedtools intersect -a sample.peaks.narrowPeak -b promoters.bed -u | wc -l)
echo "Promoter DHSs: $promoter ($(echo "scale=1; $promoter*100/$total" | bc)%)"

# Overlap with known enhancers
enhancer=$(bedtools intersect -a sample.peaks.narrowPeak -b enhancers.bed -u | wc -l)
echo "Enhancer DHSs: $enhancer ($(echo "scale=1; $enhancer*100/$total" | bc)%)"
```
