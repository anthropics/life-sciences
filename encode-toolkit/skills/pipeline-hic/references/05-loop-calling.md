# Loop Calling with HiCCUPS

HiCCUPS (Hi-C Computational Unbiased Peak Search) identifies chromatin loops
from contact matrices. It is the ENCODE-standard loop caller, part of
Juicer tools.

## HiCCUPS Loop Calling

```bash
java -Xmx13g -jar juicer_tools.jar hiccups \
    --cpu \
    --threads 4 \
    -k KR \
    -r 5000,10000,25000 \
    -f 0.1,0.1,0.1 \
    -p 4,2,1 \
    -i 7,5,3 \
    -d 20000,20000,50000 \
    sample.hic \
    loops_output/
```

This is the command `main.nf` runs with the default `--hiccups_resolutions`
(`5000,10000,25000`). The heap is 85% of the task's memory allocation, so it is
13 GB on the first attempt of the 16 GB request in `nextflow.config` and grows
with each retry, rather than being a fixed number. `--cpu` is
required with the container image, which has no CUDA runtime; drop it only with
`--hiccups_gpu` on a host with an NVIDIA GPU. CPU mode restricts the search to a
band along the diagonal (8 Mb by default), so very long-range loops are not
reported.

### Key Parameters

`-r`, `-f`, `-p`, `-i` and `-d` take one value per resolution, in the order given to
`--hiccups_resolutions`; the peak widths, window widths and merge radii are Juicer's
published defaults for 5 kb, 10 kb and 25 kb. With `--hiccups_resolutions 10000` the
workflow runs `-r 10000 -f 0.1 -p 2 -i 5 -d 20000`.

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `--cpu` | flag | CPU mode; needed without a CUDA runtime, searches 8 Mb from the diagonal |
| `-k` | KR | Normalization vector to read from the .hic file |
| `-r` | 5000,10000,25000 | Resolutions to search for loops (`--hiccups_resolutions`) |
| `-f` | 0.1,0.1,0.1 | FDR threshold per resolution |
| `-p` | 4,2,1 | Peak width (pixels) per resolution |
| `-i` | 7,5,3 | Window width (pixels) of the local background region per resolution |
| `-d` | 20000,20000,50000 | Merge radius (bp) around a loop centroid per resolution: 20 kb at 5 kb and 10 kb, 50 kb at 25 kb. juicer_tools 2.20.00 reads one value per `-r` resolution (`HiCCUPSConfiguration.extractIntegerValues(..., resolutions.length)`); its usage text still says "three values", but a list of any other length (except a single value, which is applied to every resolution) stops HiCCUPS with "Must pass N parameters" and exit code 30 |

`juicer_tools pre` must have written the `-k` vector into the .hic file:
`main.nf` builds KR, VC and VC_SQRT.

### Resolution Selection

| Resolution | Detects | Minimum Contacts |
|------------|---------|------------------|
| 5 kb | Fine-scale loops | >1 billion |
| 10 kb | Standard loops | >500 million |
| 25 kb | Large-scale loops | >100 million |

## HiCCUPS Output Format

HiCCUPS produces a BEDPE-like file with loop anchors. juicer_tools 2.20.00
writes the BEDPE core columns first, then the loop attributes in alphabetical
order, on a header line that starts with `#`:

```
#chr1  x1  x2  chr2  y1  y2  name  score  strand1  strand2  color  centroid1  centroid2  expectedBL  expectedDonut  expectedH  expectedV  fdrBL  fdrDonut  fdrH  fdrV  numCollapsed  observed  radius
```

Key columns (1-based, in that order):
- `chr1 x1 x2` (1-3) -- Upstream anchor
- `chr2 y1 y2` (4-6) -- Downstream anchor
- `centroid1`, `centroid2` (12-13) -- Centroid of the merged pixel cluster on each side
- `expectedDonut` (15) -- Expected count from the donut background model
- `fdrDonut` (19) -- FDR from the donut model (primary significance)
- `numCollapsed` (22) -- Enriched pixels merged into this call
- `observed` (23) -- Observed contact count
- `radius` (24) -- Radius of the merged cluster

The header line is the only `#` line; skip it before any arithmetic on the
file.

## Merge Loops Across Resolutions

HiCCUPS calls loops at each resolution independently and merges them itself:
it writes `merged_loops.bedpe` into the output directory alongside the
per-resolution files. There is no separate merge command (`juicer_tools` has
no `hiccups_merge` tool).

`main.nf` copies that merged file to `loops/{sample}.hiccups_loops.bedpe`,
which is the only loop file published; the per-resolution files stay in the
Nextflow work directory.

## Alternative: Mustache Loop Caller (manual, not in the image)

Mustache (Roayaei Ardakany 2020, ~165 citations) uses a scale-space
representation for loop detection. It is not installed in the container
image and the workflow does not run it:

```bash
mustache \
    -f sample.mcool \
    -r 10000 \
    -ch hg38.chrom.sizes \
    -o mustache_loops.bedpe \
    -pt 0.05 \
    -st 0.8
```

### HiCCUPS vs Mustache

| Feature | HiCCUPS | Mustache |
|---------|---------|----------|
| Background model | Donut + 3 others | Scale-space Gaussian |
| GPU support | Yes (CUDA; the image ships CPU mode only) | No |
| Speed | Faster with GPU | Moderate |
| Sensitivity | Standard | Higher (more loops) |
| ENCODE standard | Yes | Alternative |
| Concordance | ~50% overlap between callers (Wolff 2022) |

## Loop QC Metrics

### Loop Count by Resolution

HiCCUPS writes two files per resolution: `enriched_pixels_<res>.bedpe`
(pre-filter candidates) and `postprocessed_pixels_<res>.bedpe` (the final
calls that go into `merged_loops.bedpe`). Count the post-filter file:

```bash
# One file per resolution in --hiccups_resolutions (default: 5000 10000 25000)
for res in 5000 10000 25000; do
    count=$(grep -vc '^#' loops_output/postprocessed_pixels_${res}.bedpe)
    echo "Resolution ${res}: ${count} loops"
done
```

These per-resolution files are not published by the workflow; they remain in
the task's work directory. From a finished run, count the merged file instead,
excluding its header line:

```bash
grep -vc '^#' results/loops/sample.hiccups_loops.bedpe
```

Expected loop counts (human cell line, >1B contacts):
- 5 kb: 5,000-15,000 loops
- 10 kb: 3,000-10,000 loops
- 25 kb: 1,000-5,000 loops

### Loop Size Distribution

```bash
# Loop size = y1 - x1, i.e. column 5 minus column 2; skip the header line
awk '!/^#/ {print $5 - $2}' results/loops/sample.hiccups_loops.bedpe | \
    sort -n | \
    awk '{a[NR]=$1} END {
        print "Median loop size:", a[int(NR/2)];
        print "Min:", a[1];
        print "Max:", a[NR]
    }'
```

Typical loop sizes:
- Median: 200-400 kb
- Range: 50 kb to 5 Mb
- Loops <50 kb may be artifacts at lower resolutions

### CTCF Enrichment at Anchors (manual, not run by this workflow)

True loops are enriched for CTCF binding at anchors. bedtools is not in the
container image, so run this outside the workflow:

```bash
bedtools intersect \
    -a loop_anchors.bed \
    -b CTCF_peaks.bed \
    -u | wc -l
```

Expect >60% of loop anchors to overlap CTCF peaks for convergent CTCF loops.
