---
name: pipeline-cutandrun
description: "Execute CUT&RUN processing pipeline from FASTQ to peaks and signal tracks. Child of pipeline-guide. Provides Nextflow execution with Docker and cloud deployment. Use when processing CUT&RUN or CUT&Tag data, an alternative to ChIP-seq with lower background. Trigger on: CUT&RUN pipeline, CUT&Tag, SEACR, Henikoff, targeted chromatin, pA-MNase, process CUT&RUN."
---

# ENCODE CUT&RUN Pipeline: FASTQ to Peaks and Signal Tracks

## When to Use

- User wants to run a CUT&RUN or CUT&Tag processing pipeline from FASTQ to peaks
- User asks about "CUT&RUN pipeline", "CUT&Tag", "SEACR", "spike-in normalization", or "targeted chromatin"
- User needs to process CUT&RUN/CUT&Tag data with spike-in calibration and SEACR peak calling
- Example queries: "process my CUT&RUN FASTQs", "run SEACR on CUT&Tag data", "normalize CUT&RUN with spike-in controls"

Execute the CUT&RUN/CUT&Tag processing pipeline for targeted chromatin profiling,
producing peak calls with SEACR and spike-in normalized signal tracks.

## Pipeline Overview

```
FASTQ
  |-> FastQC (raw reads)
  +-> Trim Galore -> Bowtie2 (genome) -> {sample}.sorted.bam
        |
        |-> unmapped read pairs -> Bowtie2 (spike-in) -> counts -> scale_factors.txt
        |                                                                |
        +-> filter (MAPQ 10, proper pairs) -> Picard MarkDuplicates      |
            (removed) -> blacklist filter -> {sample}.filtered.bam       |
                 |-> fragment BED -> fragment bedGraph -> SEACR peaks    |
                 |-> MACS2 peaks (with --peak_caller macs2|both)         |
                 |-> FRiP vs every peak set -> {sample}.frip_mqc.tsv     |
                 +-> bamCoverage -> {sample}.normalized.bw <-- factor ---+
```

### Not run by this workflow

- **Peak-level filtering**: `--blacklist` is applied to the BAM only. Peak
  files are never filtered afterwards, and there is no separate suspect-list
  input. Pass a pre-merged blacklist + suspect-list BED as `--blacklist`, or
  filter the peak files yourself.
- **Spike-in scaling of the SEACR input**: only the bigWig is scaled
  (`bamCoverage --scaleFactor`). The fragment bedGraph given to SEACR is
  unscaled.

### ENCODE Repository

- ENCODE does not publish an official CUT&RUN pipeline. This workflow follows the published
  CUT&RUN/CUT&Tag processing protocol (Bowtie2, fragment bedGraphs, SEACR) and applies
  ENCODE conventions for filtering, blacklisting, and QC.
- **Container**: built from `scripts/Dockerfile` in this skill (`docker build -t encode-toolkit/pipeline-cutandrun:1.0.0 scripts/`); override with `--container`
- **This skill**: Nextflow DSL2 reimplementation for portability

## Core Tools and Versions

Versions are those installed by `scripts/Dockerfile`, which is what the
workflow runs.

| Tool | Version | Purpose | Citation |
|------|---------|---------|----------|
| Bowtie2 | 2.5.4 | Alignment (genome + spike-in) | Langmead & Salzberg 2012 |
| SEACR | 1.3 | Peak calling (CUT&RUN-specific) | Meers et al. 2019 |
| MACS2 | 2.2.9.1 | Alternative peak caller | Zhang et al. 2008 |
| Picard | 3.1.1 | Duplicate marking and removal | Broad Institute |
| samtools | 1.19 | BAM operations | Li et al. 2009 |
| bedtools | 2.31.0 | Genomic arithmetic | Quinlan & Hall 2010 |
| deepTools | 3.5.5 | Signal track generation | Ramirez et al. 2016 |
| Trim Galore | 0.6.10 | Adapter trimming | Krueger (Babraham) |
| FastQC | 0.12.1 | Read quality | Andrews (Babraham) |
| MultiQC | 1.21 | Aggregated QC | Ewels et al. 2016 |

The conda alternative (`cutandrun-env.yml`) pins the same version of every tool
in the table above, but installs no SEACR (only `r-base`): SEACR is not a conda
package, so on that route `SEACR_1.3.sh` and `SEACR_1.3.R` must be fetched
separately from the SEACR repository.

## Key Literature

1. **Skene & Henikoff 2017** - "An efficient targeted nuclease strategy for
   high-resolution mapping of DNA binding sites" (eLife, ~1,500 citations)
   DOI: 10.7554/eLife.21856

2. **Meers et al. 2019** - "Peak calling by Sparse Enrichment Analysis for
   CUT&RUN chromatin profiling" (Epigenetics & Chromatin, ~800 citations)
   DOI: 10.1186/s13072-019-0287-4

3. **Kaya-Okur et al. 2019** - "CUT&Tag for efficient epigenomic profiling
   of small samples and single cells" (Nature Communications, ~1,200 citations)
   DOI: 10.1038/s41467-019-09982-5

4. **Nordin et al. 2023** - "The CUT&RUN suspect list of problematic regions"
   (Genome Biology)
   DOI: 10.1186/s13059-023-02960-3

5. **Amemiya et al. 2019** - "The ENCODE Blacklist" (Scientific Reports, ~1,372 citations)
   DOI: 10.1038/s41598-019-45839-z

## Execution

### Quick Start (Local)

```bash
nextflow run scripts/main.nf \
    -profile local \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --bowtie2_index '/ref/bowtie2_index/genome' \
    --spikein_index '/ref/bowtie2_ecoli/ecoli' \
    --chrom_sizes '/ref/hg38.chrom.sizes' \
    --blacklist '/ref/hg38-blacklist.v2.bed' \
    --outdir results/ \
    -resume
```

### SLURM HPC

```bash
nextflow run scripts/main.nf \
    -profile slurm \
    --container /path/to/pipeline-cutandrun.sif \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --bowtie2_index '/ref/bowtie2_index/genome' \
    --spikein_index '/ref/bowtie2_ecoli/ecoli' \
    --chrom_sizes '/ref/hg38.chrom.sizes' \
    --blacklist '/ref/hg38-blacklist.v2.bed' \
    --outdir results/ \
    -resume
```

### Cloud (GCP / AWS)

```bash
# Google Cloud Batch
nextflow run scripts/main.nf -profile gcp \
    --container us-docker.pkg.dev/<project>/<repo>/pipeline-cutandrun:1.0.0 \
    --gcp_project <project> \
    --gcp_workdir gs://<bucket>/work \
    --reads 'gs://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --bowtie2_index gs://<bucket>/ref/bowtie2_index/genome \
    --spikein_index gs://<bucket>/ref/bowtie2_ecoli/ecoli \
    --chrom_sizes gs://<bucket>/ref/hg38.chrom.sizes \
    --blacklist gs://<bucket>/ref/hg38-blacklist.v2.bed \
    --outdir gs://<bucket>/results

# AWS Batch
nextflow run scripts/main.nf -profile aws \
    --container <account>.dkr.ecr.<region>.amazonaws.com/pipeline-cutandrun:1.0.0 \
    --aws_queue <job-queue> \
    --aws_workdir s3://<bucket>/work \
    --reads 's3://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --bowtie2_index s3://<bucket>/ref/bowtie2_index/genome \
    --spikein_index s3://<bucket>/ref/bowtie2_ecoli/ecoli \
    --chrom_sizes s3://<bucket>/ref/hg38.chrom.sizes \
    --blacklist s3://<bucket>/ref/hg38-blacklist.v2.bed \
    --outdir s3://<bucket>/results
```

`--outdir` only sets where results are published; Google Batch and AWS Batch
stage every task through the work directory, and the workflow stops with an
error if it or the project/queue is missing.

## Resource Requirements

| Step | CPUs | RAM | Time (per sample) |
|------|------|-----|-------------------|
| Bowtie2 align (genome) | 8 | 8 GB | 30-60 min |
| Bowtie2 align (spike-in) | 4 | 4 GB | 10-20 min |
| Filter/dedup | 4 | 8 GB | 15-30 min |
| SEACR peaks | 2 | 4 GB | 10-20 min |
| Signal tracks | 4 | 8 GB | 15-30 min |
| **Total** | **8** | **8 GB** | **1.5-3 hours** |

The RAM column is each step's first-attempt request. These processes ask for
that much memory per attempt, so a task killed for exceeding it is retried with
more (at most two retries, capped by `--max_memory`). Failures with any other
exit status stop the run.

## Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | required | Glob pattern to paired FASTQ files |
| `--bowtie2_index` | required | Bowtie2 genome index prefix (every file starting with this prefix is staged) |
| `--spikein_index` | `null` | Bowtie2 spike-in index prefix (E. coli by convention). When given, signal tracks are spike-in calibrated |
| `--chrom_sizes` | required | Chromosome sizes file |
| `--blacklist` | required | Blacklist BED applied to the BAM. Pass a pre-merged blacklist + CUT&RUN suspect list here if you want both |
| `--outdir` | `./results` | Output directory |
| `--seacr_mode` | `stringent` | SEACR mode: `stringent`, `relaxed`, or `both` |
| `--seacr_norm` | `norm` | SEACR normalization to the control: `norm` or `non`. Only used with `--control`; without a control the workflow always passes `non` |
| `--seacr_threshold` | `0.01` | Top fraction of signal kept by SEACR when no `--control` is given |
| `--control` | `null` | IgG control BAM, already filtered and deduplicated. Converted to a fragment bedGraph for SEACR and passed as `-c` to MACS2 |
| `--macs2_gsize` | `hs` | MACS2 effective genome size (`hs`, `mm`, or a number) |
| `--peak_caller` | `seacr` | Peak caller: `seacr`, `macs2`, or `both` |
| `--skip_spikein` | `false` | Skip spike-in calibration; signal tracks are then RPKM-normalized |

### Infrastructure parameters (`nextflow.config`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--container` | `encode-toolkit/pipeline-cutandrun:1.0.0` | Image built from `scripts/Dockerfile`. Pass a registry image for `gcp`/`aws`, or a `.sif` file for `slurm` |
| `--max_cpus`, `--max_memory`, `--max_time` | `16`, `16.GB`, `12.h` | Upper bounds applied to every process |
| `--slurm_queue`, `--slurm_account` | `normal`, none | SLURM partition and account |
| `--gcp_project`, `--gcp_workdir` | none (both required for `-profile gcp`) | Google Cloud project and `gs://` work directory |
| `--gcp_location`, `--gcp_disk` | `us-central1`, `200.GB` | Google Batch region and per-task disk |
| `--aws_queue`, `--aws_workdir` | none (both required for `-profile aws`) | AWS Batch job queue and `s3://` work directory |
| `--aws_region`, `--aws_cli_path` | `us-east-1`, `/home/ec2-user/miniconda/bin/aws` | AWS region, and the AWS CLI path inside the Batch AMI |

## Output Files

```
results/
  fastqc/                             # Raw read quality
  trim_galore/                        # Trimmed reads, trimming reports,
                                      #   and FastQC of the trimmed reads
  alignment/
    {sample}.filtered.bam             # Quality-filtered, deduplicated, blacklist-filtered
    {sample}.filtered.bam.bai
    {sample}.dup_metrics.txt          # Picard MarkDuplicates metrics
    {sample}.flagstat.txt             # samtools flagstat on the filtered BAM
  spikein/                            # Only with --spikein_index
    {sample}.spikein_counts.txt       # sample, spike-in read count
    scale_factors.txt                 # One file for the run: sample, spikein_count, scale_factor
  peaks/
    {sample}.seacr.stringent.bed      # SEACR stringent peaks
    {sample}.seacr.relaxed.bed        # With --seacr_mode relaxed or both
    {sample}.macs2_peaks.narrowPeak   # With --peak_caller macs2 or both
  signal/
    {sample}.normalized.bw            # Spike-in scaled, or RPKM without spike-in
    {sample}.fragments.bed            # Fragment BED (same chromosome, <1 kb)
  qc/
    {sample}.fragment_sizes.txt
    {sample}.frip_mqc.tsv             # FRiP, one row per peak set called for the sample
  multiqc/
    multiqc_report.html
  pipeline_info/
    timeline.html
    report.html
    trace.txt
```

The fragment bedGraph that SEACR consumes is an intermediate and is not
published; the published `signal/{sample}.fragments.bed` is the BED it is
built from.

## QC Thresholds

This is the only QC threshold table for this skill; the reference files point
back to it.

| Metric | Pass | Warning | Fail | Computed from |
|--------|------|---------|------|---------------|
| Mapping rate (genome) | >80% | 60-80% | <60% | Bowtie2 log (in `multiqc_report.html`) |
| Spike-in reads | 1-10% of total | 0.1-1% or 10-30% | <0.1% or >30% | `spikein/{sample}.spikein_counts.txt` |
| Duplication rate | <20% | 20-40% | >40% | `alignment/{sample}.dup_metrics.txt` |
| FRiP (peaks) | >10% | 5-10% | <5% | `qc/{sample}.frip_mqc.tsv` (also a MultiQC table) |
| Peak count | >5,000 | 1,000-5,000 | <1,000 | `peaks/{sample}.seacr.*.bed` |
| Fragment size | Nucleosomal pattern | Irregular | No pattern | `qc/{sample}.fragment_sizes.txt` |

### Fragment Size Distribution

CUT&RUN produces a characteristic nucleosomal ladder:
- **<120 bp**: Sub-nucleosomal (TF binding)
- **~150 bp**: Mononucleosomal (histone marks)
- **~300 bp**: Dinucleosomal
- Absence of nucleosomal pattern suggests protocol issues

## Spike-in Normalization

Spike-in normalization is CRITICAL for CUT&RUN quantitative comparison.

### How It Works

1. E. coli DNA is carried over from pA-MNase/pA-Tn5 production
2. Each sample has a different amount of spike-in reads
3. Samples with more target cleavage have fewer spike-in reads (proportionally)
4. Scale factor = smallest non-zero spike-in count across samples / this sample's count

### Scale Factor Calculation

With three samples whose spike-in counts are 200,000, 400,000 and 100,000, the
minimum is 100,000:

```
Sample A: 200,000 spike-in reads -> scale = 100,000 / 200,000 = 0.5
Sample B: 400,000 spike-in reads -> scale = 100,000 / 400,000 = 0.25
Sample C: 100,000 spike-in reads -> scale = 100,000 / 100,000 = 1.0 (minimum)
```

Higher spike-in counts = less target enrichment = lower scale factor.

All samples are written to one `spikein/scale_factors.txt` (columns: sample,
spike-in count, scale factor). A sample with no spike-in reads cannot be
calibrated and is left unscaled (factor 1). The factor is applied only to the
bigWig via `bamCoverage --scaleFactor`; the fragment bedGraph SEACR reads is
unscaled.

## SEACR vs MACS2

| Feature | SEACR | MACS2 |
|---------|-------|-------|
| Designed for | CUT&RUN/CUT&Tag | ChIP-seq |
| Background model | Sparse enrichment | Dynamic Poisson |
| Control required | Optional (IgG) | Recommended |
| Low background | Handles well | May overcall |
| Stringent mode | Very conservative | Via q-value |
| ENCODE recommendation | Primary for CUT&RUN | Alternative |

SEACR is specifically designed for the sparse, low-background signal
profile of CUT&RUN data. MACS2 may overcall peaks due to the low background.

## Critical Pitfalls

### Spike-in Calibration is CRITICAL
Without spike-in normalization, quantitative comparisons between samples are
unreliable. The amount of pA-MNase (or pA-Tn5) varies between experiments,
and spike-in reads provide the internal calibration standard. Without
`--spikein_index` (or with `--skip_spikein`) the bigWigs fall back to RPKM,
which is not quantitatively comparable across samples.

### IgG Control vs No-Antibody Control
- **IgG control**: Non-specific antibody, captures background binding
- **No-antibody**: No antibody, captures MNase accessibility background
- IgG is preferred but not always available
- SEACR can work without a control: it then uses `--seacr_threshold` (default
  0.01, the top 1% of signal) and, as SEACR v1.3 requires with a numeric
  threshold, the `non` normalization mode

### SEACR Stringent vs Relaxed Mode
- **Stringent**: Returns only the most enriched peaks (fewer, higher confidence)
- **Relaxed**: Returns a broader set including weaker peaks
- For initial analysis, use stringent mode (the default)
- For comprehensive catalogs, use `--seacr_mode both` and filter downstream

### CUT&RUN Suspect List (Nordin 2023)
The workflow applies `--blacklist` to the BAM only; it never filters the peak
files and takes no separate suspect list. To use the CUT&RUN suspect list
(Nordin et al. 2023), which identifies regions with artifactual signal
specific to CUT&RUN/CUT&Tag protocols, either pass a merged BED as
`--blacklist` or filter the peaks afterwards yourself:

```bash
# Download suspect list
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/CUTandRUN.suspectlist.hg38.bed.gz

# Option 1: merge once and pass as --blacklist (filters the BAM)
zcat CUTandRUN.suspectlist.hg38.bed.gz | cat hg38-blacklist.v2.bed - \
    | sort -k1,1 -k2,2n | bedtools merge > combined_blacklist.bed

# Option 2: filter the published peaks afterwards (manual)
bedtools intersect \
    -a results/peaks/sample.seacr.stringent.bed \
    -b combined_blacklist.bed \
    -v \
    > sample_peaks_filtered.bed
```

### CUT&RUN vs CUT&Tag
Both protocols are supported by this pipeline. Differences:
- **CUT&RUN**: Uses pA-MNase, E. coli spike-in from MNase production
- **CUT&Tag**: Uses pA-Tn5, E. coli spike-in from Tn5 production
- CUT&Tag has higher background from Tn5 insertion preference
- CUT&Tag may work better for histone marks; CUT&RUN for TFs

## Provenance Integration

After pipeline completion, log all outputs:

```python
encode_log_derived_file(
    file_path="/results/peaks/sample1.seacr.stringent.bed",
    source_accessions=["ENCSR...", "ENCFF..."],
    description="CUT&RUN peaks from ENCODE CUT&RUN pipeline",
    file_type="CUT&RUN_peaks",
    tool_used="Bowtie2 2.5.4 + SEACR 1.3",
    parameters="stringent mode, threshold 0.01 non, BAM blacklist-filtered (peaks unfiltered)"
)
```

## Reference Files

Detailed step-by-step documentation is provided in the `references/` directory:

1. `01-qc-trimming.md` -- Read QC and adapter trimming for CUT&RUN
2. `02-bowtie2-alignment.md` -- Bowtie2 alignment to genome and spike-in
3. `03-filtering-spikein.md` -- Filtering, dedup, and spike-in normalization
4. `04-seacr-peaks.md` -- SEACR peak calling and MACS2 alternative
5. `05-qc-metrics.md` -- Fragment sizes, FRiP, spike-in QC

## Walkthrough: Processing ENCODE CUT&RUN from FASTQ to Peaks

**Goal**: Process CUT&RUN/CUT&Tag FASTQ files through the ENCODE-compatible pipeline to generate peak calls with spike-in normalization.
**Context**: CUT&RUN uses targeted MNase digestion (lower background than ChIP-seq) but requires different peak calling (SEACR instead of MACS2) and spike-in normalization for quantitative comparisons.

### Step 1: Find CUT&RUN experiment

```
encode_search_experiments(assay_title="CUT&RUN", organism="Homo sapiens")
```

Expected output:
```json
{
  "results": [
    {"accession": "ENCSR900CUR", "assay_title": "CUT&RUN", "target": "H3K27me3", "biosample_summary": "K562", "assembly": ["GRCh38"], "status": "released"}
  ],
  "total": 35,
  "limit": 25,
  "offset": 0,
  "has_more": true,
  "next_offset": 25
}
```

### Step 2: List FASTQ files

```
encode_list_files(experiment_accession="ENCSR900CUR", file_format="fastq")
```

Expected output (a JSON array of file records; fields abridged):
```json
[
  {"accession": "ENCFF900CR1", "file_format": "fastq", "output_type": "reads", "biological_replicates": [1], "file_size": 839252000, "file_size_human": "800.4 MB", "status": "released"},
  {"accession": "ENCFF901CR2", "file_format": "fastq", "output_type": "reads", "biological_replicates": [1], "file_size": 891394000, "file_size_human": "850.1 MB", "status": "released"}
]
```

**Interpretation**: CUT&RUN yields smaller files than ChIP-seq (~800MB vs ~2.5GB) due to lower background.

### Step 3: Name the files so a read-pair glob can find them

ENCODE FASTQs are named by accession, so the two mates of a pair share no
prefix, and the workflow matches file pairs with a `{1,2}` glob. Which mate a
file is comes from its page on encodeproject.org (`paired_end` 1 or 2, and
`paired_with` naming the other accession), not from any tool here. Link the
files into the shape the glob expects:

```bash
mkdir -p fastq
ln -s "$PWD/ENCFF900CR1.fastq.gz" fastq/ENCSR900CUR_R1.fastq.gz
ln -s "$PWD/ENCFF901CR2.fastq.gz" fastq/ENCSR900CUR_R2.fastq.gz
```

### Step 4: Run the CUT&RUN pipeline

```bash
nextflow run scripts/main.nf \
  -profile local \
  --reads 'fastq/ENCSR900CUR_R{1,2}.fastq.gz' \
  --bowtie2_index '/ref/bowtie2_index/genome' \
  --spikein_index '/ref/bowtie2_ecoli/ecoli' \
  --chrom_sizes '/ref/hg38.chrom.sizes' \
  --blacklist '/ref/hg38-blacklist.v2.bed' \
  --peak_caller seacr \
  --outdir results/ \
  -resume
```

Key pipeline steps:
1. FastQC on the raw reads, then adapter trimming (Trim Galore, `--nextera`)
2. Bowtie2 alignment (`--very-sensitive --no-mixed --no-discordant --dovetail -I 10 -X 700`)
3. Spike-in alignment of the read pairs that did not map to the genome
4. Scale factor per sample (minimum spike-in count / sample count)
5. Filter (MAPQ 10, proper pairs), remove duplicates with Picard, remove blacklist regions
6. SEACR peak calling from the fragment bedGraph (stringent by default)
7. Signal bigWig with `bamCoverage`, scaled by the spike-in factor
8. FRiP of the filtered BAM against every peak set called for the sample,
   written to `qc/{sample}.frip_mqc.tsv`

### Step 5: Validate output quality

Use the QC threshold table above with `alignment/{sample}.dup_metrics.txt`,
`spikein/{sample}.spikein_counts.txt`, `qc/{sample}.fragment_sizes.txt` and
`qc/{sample}.frip_mqc.tsv`.

**Key difference from ChIP-seq**: CUT&RUN has inherently lower background, so peak callers like MACS2 overfit. Use SEACR (Meers et al. 2019) instead.

### Step 6: Compare with ChIP-seq for the same target

```
encode_search_experiments(assay_title="Histone ChIP-seq", biosample_term_name="K562", target="H3K27me3", organism="Homo sapiens")
```

**Interpretation**: CUT&RUN typically identifies fewer but higher-confidence peaks than ChIP-seq. Concordant peaks between both methods are the highest confidence.

### Integration with downstream skills
- SEACR peaks feed into -> **histone-aggregation** for cross-experiment comparison
- Spike-in normalized signals feed into -> **visualization-workflow**
- Peak regions feed into -> **regulatory-elements** for chromatin state classification
- QC uses different thresholds than ChIP-seq -> **quality-assessment** (see suspect list)
- Pipeline provenance logged by -> **data-provenance**

## Code Examples

### 1. Survey CUT&RUN/CUT&Tag availability

```
encode_get_facets(assay_title="CUT&RUN", organism="Homo sapiens")
```

Expected output:
```json
{
  "target.label": [
    {"term": "H3K27me3", "count": 15},
    {"term": "H3K4me3", "count": 12},
    {"term": "H3K27ac", "count": 8},
    {"term": "CTCF", "count": 5}
  ]
}
```

### 2. Find matching ChIP-seq for comparison

```
encode_search_experiments(assay_title="Histone ChIP-seq", biosample_term_name="K562", target="H3K27me3", organism="Homo sapiens")
```

Expected output:
```json
{
  "results": [
    {"accession": "ENCSR000CHI", "assay_title": "Histone ChIP-seq", "target": "H3K27me3", "biosample_summary": "K562", "assembly": ["GRCh38"]}
  ],
  "total": 5,
  "limit": 25,
  "offset": 0,
  "has_more": false,
  "next_offset": null
}
```

### 3. Track CUT&RUN experiments

```
encode_track_experiment(accession="ENCSR900CUR", notes="K562 H3K27me3 CUT&RUN - SEACR peaks for comparison with ChIP-seq")
```

Expected output:
```json
{
  "tracking": {
    "accession": "ENCSR900CUR",
    "action": "tracked"
  },
  "publications_found": 0,
  "publications": [],
  "pipelines_found": 0,
  "pipelines": []
}
```

## Integration

| This skill produces... | Feed into... | Purpose |
|---|---|---|
| SEACR peaks | **histone-aggregation** | Cross-experiment comparison (note: different caller than ChIP-seq) |
| Spike-in normalized signal | **visualization-workflow** | Quantitatively comparable browser tracks |
| Peak regions | **regulatory-elements** | Chromatin state classification |
| CUT&RUN-specific QC | **quality-assessment** | Validate with CUT&RUN-appropriate thresholds |
| Peak coordinates | **motif-analysis** | TF motif discovery at CUT&RUN peaks |
| Pipeline parameters | **data-provenance** | Record SEACR/spike-in normalization details |
| Peak files | **variant-annotation** | Identify variants in CUT&RUN peaks |
| Comparison with ChIP-seq | **compare-biosamples** | Cross-assay concordance analysis |

## Related Skills

- `pipeline-guide` -- Parent skill with compute resource assessment and cloud setup
- `histone-aggregation` -- Aggregate histone mark data across samples
- `quality-assessment` -- Evaluate pipeline output quality metrics
- `data-provenance` -- Track all pipeline inputs, outputs, and parameters
- `download-encode` -- Download ENCODE CUT&RUN FASTQ files for pipeline input
- `publication-trust` -- Verify literature claims backing analytical decisions

## Presenting Results

When reporting CUT&RUN pipeline results:

- **SEACR peak counts**: Report peak counts for each SEACR mode that was run (default: stringent only; both with `--seacr_mode both`). If MACS2 was also run, include those counts for comparison
- **Spike-in normalization factor**: Report the scale factor and spike-in count per sample from `spikein/scale_factors.txt` and the spike-in read fraction (ideal 1-10% of total reads). Explain that higher spike-in counts indicate less target enrichment
- **FRiP**: Report it from `qc/{sample}.frip_mqc.tsv`, which has one row per peak set called for the sample (SEACR stringent and/or relaxed, and/or MACS2), and judge each against the QC table (>10% pass, 5-10% warning, <5% fail). The value is the fraction of the filtered BAM's alignments that overlap a peak, so mates of a pair count separately; the `Peak set` column names the peak file each row refers to
- **Signal track paths**: Provide paths to the `signal/{sample}.normalized.bw` files (spike-in scaled, or RPKM if spike-in was skipped) for genome browser visualization
- **Fragment size distribution**: From `qc/{sample}.fragment_sizes.txt`, confirm the expected nucleosomal ladder pattern and note the dominant fragment class (sub-nucleosomal for TFs, mononucleosomal for histone marks)
- **Key QC metrics**: Present mapping rate (>80%), duplication rate (<20%), and spike-in calibration status in a summary table
- **Blacklist filtering**: State that `--blacklist` was applied to the BAM and that the peak files are unfiltered; note separately whether a suspect list was merged into `--blacklist` or applied to the peaks manually
- **Next steps**: Suggest `peak-annotation` for gene association of peaks, or `visualization-workflow` for genome browser session generation

## For the request: "$ARGUMENTS"
