---
name: pipeline-dnaseseq
description: "Execute ENCODE DNase-seq pipeline from FASTQ to hotspots and footprints. Child of pipeline-guide. Provides Nextflow execution with Docker and cloud deployment. Use when processing DNase-seq data, calling DNase hypersensitive sites, performing footprinting analysis. Trigger on: DNase-seq pipeline, DNase hypersensitive, DHS, Hotspot2, footprinting, DNase I, chromatin accessibility DNase."
---

# ENCODE DNase-seq Pipeline: FASTQ to Hotspots and Footprints

## When to Use

- User wants to run a DNase-seq processing pipeline from FASTQ to hotspots and footprints
- User asks about "DNase-seq pipeline", "DNase hypersensitive sites", "Hotspot2", "footprinting", or "DHS"
- User needs to process DNase-seq data for chromatin accessibility and TF footprint analysis
- Example queries: "process my DNase-seq FASTQs", "call DNase hypersensitive sites", "run footprinting analysis on DNase-seq"

Execute the ENCODE DNase-seq pipeline for chromatin accessibility profiling,
producing DNase hypersensitive sites (DHSs) via Hotspot2 and transcription
factor footprints.

## Pipeline Overview

```
FASTQ -> Trim -> BWA-MEM align -> Filter/dedup -> Hotspot2 -> DHS peaks
                                       |                        |
                                    Signal track         Footprinting (HINT)
```

### ENCODE Repository

- **GitHub**: `ENCODE-DCC/dnase-seq-pipeline`
- **Container**: built from `scripts/Dockerfile` in this skill (`docker build -t encode-toolkit/pipeline-dnaseseq:1.0.0 scripts/`); override with `--container`
- **WDL**: Available for Cromwell execution
- **This skill**: Nextflow DSL2 reimplementation for portability

## Core Tools and Versions

| Tool | Version | Purpose | Citation |
|------|---------|---------|----------|
| BWA-MEM | 0.7.18 | Alignment | Li & Durbin 2009 |
| samtools | 1.19 | BAM operations | Li et al. 2009 |
| Picard | 3.1.1 | Duplicate marking | Broad Institute |
| Hotspot2 | 2.1.2 | DHS calling (ENCODE standard) | John et al. 2011 |
| modwt | 1.0 | Wavelet smoothing used by Hotspot2 | Stam Lab |
| bedtools | 2.31.0 | Genomic arithmetic | Quinlan & Hall 2010 |
| BEDOPS | apt (Ubuntu 22.04) | `sort-bed` and `unstarch` for the Hotspot2 starch archives | Neph et al. 2012 |
| HINT (RGT) | 1.0.2 | TF footprinting | Li et al. 2019 |
| FastQC | 0.12.1 | Read quality | Andrews (Babraham) |
| Trim Galore | 0.6.10 | Adapter and quality trimming | Krueger (Babraham) |
| cutadapt | 4.6 | Adapter removal backend for Trim Galore | Martin 2011 |
| MultiQC | 1.21 | Aggregated QC | Ewels et al. 2016 |

## Key Literature

1. **John et al. 2011** - "Chromatin accessibility pre-determines glucocorticoid
   receptor binding patterns" (Nature Genetics, ~600 citations)
   DOI: 10.1038/ng.759

2. **Thurman et al. 2012** - "The accessible chromatin landscape of the human
   genome" (Nature, ~3,000 citations)
   DOI: 10.1038/nature11232

3. **Vierstra et al. 2020** - "Global reference mapping of human transcription
   factor footprints" (Nature, ~600 citations)
   DOI: 10.1038/s41586-020-2528-x

4. **Amemiya et al. 2019** - "The ENCODE Blacklist" (Scientific Reports, ~1,372 citations)
   DOI: 10.1038/s41598-019-45839-z

5. **Li et al. 2019** - "Identification of transcription factor binding sites using
   ATAC-seq" (Genome Biology) -- HINT-ATAC footprinting
   DOI: 10.1186/s13059-019-1642-2

## Execution

### Quick Start (Local)

```bash
nextflow run scripts/main.nf \
    -profile local \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index '/ref/bwa_index/genome.fa' \
    --chrom_sizes '/ref/hg38.chrom.sizes' \
    --hotspot_center_sites '/ref/hotspot2/hg38.center_sites.n100.starch' \
    --hotspot_mappable '/ref/hotspot2/hg38.mappable_only.bed' \
    --rgt_data '/ref/rgtdata' \
    --blacklist '/ref/hg38-blacklist.v2.bed' \
    --outdir results/ \
    -resume
```

Drop `--rgt_data` and add `--skip_footprint` to stop after hotspot calling.

### SLURM HPC

The `slurm` profile runs through Singularity, which cannot resolve the default
Docker image name, so pass the converted `.sif` with `--container`:

```bash
singularity build pipeline-dnaseseq.sif docker-daemon://encode-toolkit/pipeline-dnaseseq:1.0.0

nextflow run scripts/main.nf \
    -profile slurm \
    --container /path/to/pipeline-dnaseseq.sif \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index '/ref/bwa_index/genome.fa' \
    --chrom_sizes '/ref/hg38.chrom.sizes' \
    --hotspot_center_sites '/ref/hotspot2/hg38.center_sites.n100.starch' \
    --hotspot_mappable '/ref/hotspot2/hg38.mappable_only.bed' \
    --rgt_data '/ref/rgtdata' \
    --blacklist '/ref/hg38-blacklist.v2.bed' \
    --outdir results/ \
    -resume
```

### Cloud (GCP / AWS)

```bash
# Google Cloud Batch
nextflow run scripts/main.nf -profile gcp \
    --container us-docker.pkg.dev/<project>/<repo>/pipeline-dnaseseq:1.0.0 \
    --gcp_project <project> \
    --gcp_workdir gs://<bucket>/work \
    --reads 'gs://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index gs://<bucket>/ref/bwa_index/genome.fa \
    --chrom_sizes gs://<bucket>/ref/hg38.chrom.sizes \
    --hotspot_center_sites gs://<bucket>/ref/hotspot2/hg38.center_sites.n100.starch \
    --hotspot_mappable gs://<bucket>/ref/hotspot2/hg38.mappable_only.bed \
    --rgt_data gs://<bucket>/ref/rgtdata \
    --blacklist gs://<bucket>/ref/hg38-blacklist.v2.bed \
    --outdir gs://<bucket>/results

# AWS Batch
nextflow run scripts/main.nf -profile aws \
    --container <account>.dkr.ecr.<region>.amazonaws.com/pipeline-dnaseseq:1.0.0 \
    --aws_queue <job-queue> \
    --aws_workdir s3://<bucket>/work \
    --reads 's3://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index s3://<bucket>/ref/bwa_index/genome.fa \
    --chrom_sizes s3://<bucket>/ref/hg38.chrom.sizes \
    --hotspot_center_sites s3://<bucket>/ref/hotspot2/hg38.center_sites.n100.starch \
    --hotspot_mappable s3://<bucket>/ref/hotspot2/hg38.mappable_only.bed \
    --rgt_data s3://<bucket>/ref/rgtdata \
    --blacklist s3://<bucket>/ref/hg38-blacklist.v2.bed \
    --outdir s3://<bucket>/results
```

`--outdir` only sets where results are published; Google Batch and AWS Batch
stage every task through the work directory, and the workflow stops with an
error if it or the project/queue is missing.

## Resource Requirements

| Step | CPUs | RAM | Time (per sample) |
|------|------|-----|-------------------|
| BWA-MEM align | 8 | 16 GB | 1-2 hours |
| Filter/dedup | 4 | 8 GB | 30-60 min |
| Hotspot2 | 4 | 8 GB | 30-60 min |
| Signal generation | 2 | 4 GB | 15-30 min |
| Footprinting | 4 | 8 GB | 1-2 hours |
| **Total** | **8** | **16 GB** | **3-6 hours** |

Every process asks for `memory { N.GB * task.attempt }`, so a task killed for running out
of memory is retried with more: the second attempt gets twice the figure in the table, the
third three times it, bounded by `--max_memory` (32 GB by default). `nextflow.config`
scales the time the same way for `BWA_ALIGN`, `FILTER_DEDUP`, `HOTSPOT2` and
`FOOTPRINTING`, bounded by `--max_time`; the other processes declare no time limit. A task
is retried only for exit codes 130-145 and 104 (killed for exceeding a limit); any other
failure stops the run.

## Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | required | Glob pattern to paired FASTQ files (e.g. `'/data/*_R{1,2}.fastq.gz'`). Paired-end only |
| `--bwa_index` | required | BWA index **prefix**, i.e. the FASTA path. Every file matching `<prefix>*` (the `.fa` plus `.amb .ann .bwt .pac .sa`) is staged |
| `--chrom_sizes` | required | Two-column chromosome sizes file. Used for the bigWig track and converted to BED for Hotspot2 `-c` |
| `--hotspot_center_sites` | required | Hotspot2 center-sites archive (`.starch`), made once per genome with `extractCenterSites.sh` (Hotspot2 `-C`) |
| `--hotspot_mappable` | `null` | Mappable-regions BED that the center sites were made from (Hotspot2 `-M`; recommended) |
| `--blacklist` | required | ENCODE blacklist BED. Applied to the BAM before hotspot calling |
| `--outdir` | `./results` | Output directory |
| `--fdr` | `0.05` | Hotspot2 hotspot FDR (`-f`). Names every Hotspot2 output file. `-F` is passed as `max(--fdr, 0.05)` because it may not be stricter than `-f` |
| `--skip_footprint` | `false` | Skip footprinting analysis |
| `--organism` | `hg38` | Genome name registered in the RGT data directory, used by HINT footprinting |
| `--rgt_data` | required unless `--skip_footprint` | RGT data directory with the genome for `--organism` set up (see below) |

There is no `--genome`, `--single_end`, or `--fastq_r1/--fastq_r2` parameter.

### Infrastructure parameters (`nextflow.config`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--container` | `encode-toolkit/pipeline-dnaseseq:1.0.0` | Image built from `scripts/Dockerfile`. Pass a registry image for `gcp`/`aws`, or a `.sif` file for `slurm` |
| `--max_cpus`, `--max_memory`, `--max_time` | `16`, `32.GB`, `24.h` | Upper bounds applied to every process |
| `--slurm_queue`, `--slurm_account` | `normal`, none | SLURM partition and account |
| `--gcp_project`, `--gcp_workdir` | none (both required for `-profile gcp`) | Google Cloud project and `gs://` work directory |
| `--gcp_location`, `--gcp_disk` | `us-central1`, `200.GB` | Google Batch region and per-task disk |
| `--aws_queue`, `--aws_workdir` | none (both required for `-profile aws`) | AWS Batch job queue and `s3://` work directory |
| `--aws_region`, `--aws_cli_path` | `us-east-1`, `/home/ec2-user/miniconda/bin/aws` | AWS region, and the AWS CLI path inside the Batch AMI |

## Output Files

```
results/
  fastqc/                          # FastQC on the raw reads
  trim_galore/
    {sample}_R1_val_1.fq.gz        # Trimmed reads
    {sample}_R2_val_2.fq.gz
    *_trimming_report.txt          # Trim Galore reports
    *_fastqc.{html,zip}            # FastQC on the trimmed reads
  alignment/
    {sample}.filtered.bam          # Filtered, deduplicated, blacklist-free BAM
    {sample}.filtered.bam.bai
    {sample}.flagstat.txt          # samtools flagstat on the filtered BAM
    {sample}.dup_metrics.txt       # Picard MarkDuplicates metrics
  hotspots/
    {sample}.hotspots.fdr0.05.bed  # DHS hotspots (primary output; unstarched)
    {sample}.peaks.narrowPeak      # Peaks within hotspots (unstarched)
    {sample}.allcalls.bed          # All site calls before FDR filtering (unstarched)
    {sample}.SPOT.txt              # SPOT score
    {sample}.density.bw            # RPM fragment-coverage signal track (bigWig)
  footprints/
    {sample}.footprints.bed        # TF footprints (omitted with --skip_footprint)
  qc/
    {sample}.insert_sizes.txt      # samtools stats output (insert sizes in the IS block)
  multiqc/
    multiqc_report.html
  pipeline_info/
    timeline.html
    report.html
    trace.txt
```

The `.bed`, `.narrowPeak` and `.SPOT.txt` files under `hotspots/` are the
`unstarch`-ed forms of the Hotspot2 `.starch` archives; the raw archives stay in
the Nextflow work directory. `{sample}.density.bw` is the bedtools RPM fragment
coverage track, not the per-base cut-count bigWig Hotspot2 writes internally.

## QC Thresholds (ENCODE Standards)

| Metric | Pass | Warning | Fail |
|--------|------|---------|------|
| SPOT score (Signal Portion of Tags) | >0.4 | 0.2-0.4 | <0.2 |
| Hotspot count | >50,000 | 20,000-50,000 | <20,000 |
| Mapping rate | >80% | 60-80% | <60% |
| Duplication rate | <30% | 30-50% | >50% |
| NRF (Non-Redundant Fraction) | >0.8 | 0.7-0.8 | <0.7 |
| PBC1 (PCR Bottleneck Coefficient 1) | >0.9 | 0.7-0.9 | <0.7 |
| Insert size peak | 50-150 bp | Variable | Abnormal |

The workflow produces everything the first four rows and the last row need:
the SPOT score (`hotspots/{sample}.SPOT.txt`), the hotspot BED to count, the
mapping and duplication rates (`alignment/{sample}.flagstat.txt` and
`{sample}.dup_metrics.txt`), and the insert-size distribution
(`qc/{sample}.insert_sizes.txt`). NRF, PBC1 and PBC2 are **not** computed;
derive them manually from the alignment BAM as shown in
`references/03-filtering.md`.

### SPOT Score

The SPOT score (Signal Portion of Tags) is the fraction of reads falling
within hotspots. It is the DNase-seq equivalent of FRiP for ChIP-seq.

Higher SPOT = more enrichment in accessible regions = better library quality.

## Hotspot2 vs MACS2

**IMPORTANT**: ENCODE uses Hotspot2 for DNase-seq, NOT MACS2.

| Feature | Hotspot2 | MACS2 |
|---------|----------|-------|
| Designed for | DNase-seq | ChIP-seq |
| Background model | Local tag density + mappability | Dynamic Poisson |
| ENCODE standard | Yes (DNase-seq) | Yes (ChIP-seq/ATAC-seq) |
| Mappability correction | Built-in | Not available |
| Output | Hotspots + peaks | Peaks only |

Hotspot2 accounts for mappability variation across the genome, which is
critical for DNase-seq because DNase I cuts accessible chromatin regardless
of whether it is uniquely mappable.

## Critical Pitfalls

### DNase-seq vs ATAC-seq
These are different assays measuring the same biology (chromatin accessibility):
- **DNase-seq**: Uses DNase I enzyme, requires more input material
- **ATAC-seq**: Uses Tn5 transposase, works on fewer cells
- Analysis pipelines differ: Hotspot2 for DNase-seq, MACS2 for ATAC-seq
- Data are largely concordant but not identical

### Fragment Size Distribution
DNase-seq produces a characteristic fragment size distribution:
- Peak at ~50-100 bp (sub-nucleosomal fragments at DHS)
- Secondary peak at ~150-200 bp (mononucleosomal fragments)
- Long tail of larger fragments
- If distribution is abnormal, check library preparation protocol

### Mappability Index
Hotspot2 needs a center-sites file, which is derived from a mappable-regions BED. Both are
read-length and genome-build specific. Create the center sites once per genome with the script
that ships with Hotspot2 (it is on the PATH inside the image):

```bash
# chrom_sizes.bed is a BED file: chromosome, 0, length
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' hg38.chrom.sizes | sort-bed - > chrom_sizes.bed
extractCenterSites.sh -c chrom_sizes.bed -M hg38.mappable_only.bed -o hg38.center_sites.n100.starch
```

Pass the same mappable-regions BED to the workflow as `--hotspot_mappable` that
was used to build the center sites.

Mappable-regions files:
- hg38 / 36 bp: Use ENCODE-provided index
- hg38 / 76 bp: Use ENCODE-provided index
- hg38 / 150 bp: May need to generate custom index
- Wrong mappability index = incorrect peak calls

The workflow ends at footprint calling; motif matching against JASPAR is a separate
downstream step (see the `jaspar-motifs` skill).

### RGT Data Directory (footprinting)
HINT reads genome sequence and annotation from an RGT data directory, which is several GB and
is not part of the image. Create it once, then pass it with `--rgt_data`:

```bash
pip install RGT==1.0.2            # creates ~/rgtdata with setupGenomicData.py
cd ~/rgtdata && python setupGenomicData.py --hg38
```

Then run with `--rgt_data ~/rgtdata`, or use `--skip_footprint` to stop after hotspot calling.

### Blacklist Filtering
`--blacklist` is required and the workflow removes blacklisted reads from the
BAM (Amemiya et al. 2019) before Hotspot2 runs, so the published peaks are
already blacklist-free. Blacklist regions produce artifactual signal in
accessibility assays.

Filter at the peak level only when the peaks came from a BAM that was not
filtered, or when applying an additional list:
```bash
bedtools intersect -a hotspots.bed -b hg38-blacklist.v2.bed -v > hotspots_filtered.bed
```

## Footprinting Analysis

Transcription factor footprinting detects bound TFs from DNase-seq signal.
This is what the workflow runs:

### HINT Footprinting (DNase-seq mode)
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

Use `--dnase-seq`, not `--atac-seq`: the two apply different cleavage-bias
models, and the wrong one silently produces wrong footprints.

### Interpretation
- Footprints are depressions in the DNase signal where a bound TF protects DNA
- Requires deep sequencing (>100M reads) for reliable footprints
- Sensitivity varies by TF: pioneer factors have shallow footprints
- Vierstra et al. 2020 provides a global reference map for comparison

## Provenance Integration

After pipeline completion, log all outputs:

```python
encode_log_derived_file(
    file_path="/results/hotspots/sample1.hotspots.fdr0.05.bed",
    source_accessions=["ENCSR...", "ENCFF..."],
    description="DNase hypersensitive sites from ENCODE DNase-seq pipeline",
    file_type="DHS_peaks",
    tool_used="BWA 0.7.18 + Hotspot2 2.1.2",
    parameters="FDR 0.05, blacklist filtered, ENCODE hg38 mappability index"
)
```

## Reference Files

Detailed step-by-step documentation is provided in the `references/` directory:

1. `01-qc-trimming.md` -- Read QC and adapter trimming
2. `02-alignment.md` -- BWA-MEM alignment for DNase-seq
3. `03-filtering.md` -- BAM filtering, deduplication, blacklist removal
4. `04-hotspot-calling.md` -- Hotspot2 DHS detection and signal generation
5. `05-footprinting.md` -- TF footprint detection with HINT

## Walkthrough: Processing ENCODE DNase-seq from FASTQ to Hypersensitive Sites

**Goal**: Process raw DNase-seq FASTQ files through the ENCODE pipeline to generate DNase I hypersensitive site (DHS) peak calls.
**Context**: DNase-seq identifies open chromatin via DNase I enzyme digestion. The pipeline uses BWA alignment and Hotspot2 for DHS identification.

### Step 1: Find DNase-seq experiment

```
encode_search_experiments(assay_title="DNase-seq", biosample_term_name="K562", organism="Homo sapiens")
```

Expected output:
```json
{
  "results": [
    {"accession": "ENCSR000DNS", "assay_title": "DNase-seq", "biosample_summary": "K562", "status": "released"}
  ],
  "total": 8,
  "limit": 25,
  "offset": 0,
  "has_more": false,
  "next_offset": null
}
```

### Step 2: List and download FASTQ files

```
encode_list_files(experiment_accession="ENCSR000DNS", file_format="fastq")
```

Expected output (a JSON array of file records; fields abridged):
```json
[
  {"accession": "ENCFF500DN1", "file_format": "fastq", "output_type": "reads", "file_size_human": "2.6 GB", "biological_replicates": [1], "status": "released"},
  {"accession": "ENCFF501DN2", "file_format": "fastq", "output_type": "reads", "file_size_human": "2.7 GB", "biological_replicates": [1], "status": "released"}
]
```

```
encode_download_files(file_accessions=["ENCFF500DN1", "ENCFF501DN2"], download_dir="/data/dnaseseq/fastq")
```

### Step 3: Name the FASTQs so a read-pair glob can find them, then run the pipeline

ENCODE names every FASTQ after its accession (`ENCFF500DN1.fastq.gz`), with no `_R1`/`_R2`
in the name, so the two files of a pair share no prefix and the `--reads` glob cannot pair
them. Link them into the shape the glob expects. Which mate an accession is comes from the
ENCODE file record on encodeproject.org, which carries `paired_end` (1 or 2) and
`paired_with`; the MCP file tools do not return those two fields:

```bash
cd /data/dnaseseq/fastq
ln -s ENCFF500DN1.fastq.gz k562_rep1_R1.fastq.gz
ln -s ENCFF501DN2.fastq.gz k562_rep1_R2.fastq.gz
```

```bash
nextflow run scripts/main.nf \
  -profile local \
  --reads '/data/dnaseseq/fastq/k562_*_R{1,2}.fastq.gz' \
  --bwa_index /ref/bwa_index/genome.fa \
  --chrom_sizes /ref/hg38.chrom.sizes \
  --hotspot_center_sites /ref/hotspot2/hg38.center_sites.n100.starch \
  --hotspot_mappable /ref/hotspot2/hg38.mappable_only.bed \
  --rgt_data /ref/rgtdata \
  --blacklist /ref/hg38-blacklist.v2.bed \
  --outdir results/ \
  -resume
```

Key pipeline steps:
1. Quality trimming
2. BWA-MEM alignment
3. Duplicate removal and blacklist filtering
4. Hotspot2 DHS calling
5. Signal track generation
6. Footprint analysis (HINT, DNase-seq mode)

### Step 4: Validate output quality

| Metric | Threshold | Purpose |
|---|---|---|
| SPOT score | > 0.4 | Signal portion of tags |
| Hotspot count | > 50,000 | Sensitivity |
| Duplicate rate | < 30% | Library complexity |

### Step 5: Compare with ATAC-seq

```
encode_search_experiments(assay_title="ATAC-seq", biosample_term_name="K562", organism="Homo sapiens")
```

**Interpretation**: DNase-seq and ATAC-seq both measure accessibility but with different biases. Compare peaks from both assays -- concordant peaks are high confidence.

### Integration with downstream skills
- DHS peaks feed into -> **accessibility-aggregation** alongside ATAC-seq peaks
- Footprint data feeds into -> **motif-analysis** for TF binding prediction
- Signal tracks feed into -> **visualization-workflow**
- Peaks integrate with -> **regulatory-elements** for cCRE classification

## Code Examples

### 1. Survey DNase-seq availability

```
encode_get_facets(assay_title="DNase-seq", organism="Homo sapiens")
```

Expected output (facet field names are the top-level keys):
```json
{
  "biosample_ontology.organ_slims": [
    {"term": "blood", "count": 45},
    {"term": "brain", "count": 30}
  ]
}
```

### 2. Check for existing DHS peaks

```
encode_list_files(experiment_accession="ENCSR000DNS", file_format="bed", output_type="peaks", assembly="GRCh38")
```

Expected output (a JSON array of file records; fields abridged):
```json
[
  {"accession": "ENCFF800DHS", "file_format": "bed", "file_type": "bed narrowPeak", "output_type": "peaks", "assembly": "GRCh38", "file_size_human": "1.5 MB"}
]
```

### 3. Track DNase-seq experiments

```
encode_track_experiment(accession="ENCSR000DNS", notes="K562 DNase-seq for accessibility comparison with ATAC-seq")
```

Expected output (the `notes` you pass are stored, not echoed back; read them with `encode_list_tracked`):
```json
{
  "tracking": {"accession": "ENCSR000DNS", "action": "tracked"},
  "publications_found": 0,
  "publications": [],
  "pipelines_found": 1,
  "pipelines": [
    {"title": "DNase-HS pipeline single-end - Version 2", "version": "2.0", "software": [{"name": "bwa", "version": "0.7.17"}], "status": "released"}
  ]
}
```

## Integration

| This skill produces... | Feed into... | Purpose |
|---|---|---|
| DHS peaks (narrowPeak) | **accessibility-aggregation** | Union merge with ATAC-seq peaks |
| TF footprints | **motif-analysis** | Validate motif predictions with footprint evidence |
| Signal tracks (bigWig) | **visualization-workflow** | Genome browser display |
| Accessible regions | **regulatory-elements** | cCRE classification |
| DHS coordinates | **variant-annotation** | Annotate variants in hypersensitive sites |
| QC metrics | **quality-assessment** | Validate SPOT score and sensitivity |
| Pipeline parameters | **data-provenance** | Record BWA/Hotspot2 versions |
| DHS peak regions | **jaspar-motifs** | Scan accessible sites for known TF motifs |

## Related Skills

- `pipeline-guide` -- Parent skill with compute resource assessment and cloud setup
- `accessibility-aggregation` -- Aggregate DHS data across samples/tissues
- `quality-assessment` -- Evaluate pipeline output quality metrics
- `data-provenance` -- Track all pipeline inputs, outputs, and parameters
- `download-encode` -- Download ENCODE DNase-seq FASTQ files for pipeline input
- `publication-trust` -- Verify literature claims backing analytical decisions

## Presenting Results

When reporting DNase-seq pipeline results:

- **Hotspot counts**: Report total Hotspot2 DHS calls at the specified FDR threshold and the number remaining after blacklist filtering
- **Signal-to-noise (SPOT score)**: Report the SPOT score prominently (>0.4 pass, 0.2-0.4 warning, <0.2 fail). This is the DNase-seq equivalent of FRiP
- **Footprint depth**: If footprinting was performed, report the number of lines in `footprints/{sample}.footprints.bed` and note the sequencing depth (>100M reads recommended for reliable footprints)
- **Key QC metrics**: Present mapping rate (>80%) and duplication rate (<30%) from `alignment/{sample}.flagstat.txt` and `alignment/{sample}.dup_metrics.txt`, and the insert size peak from the `IS` block of `qc/{sample}.insert_sizes.txt`. NRF (>0.8) and PBC1 (>0.9) are not produced by the workflow -- state that they were computed manually (`references/03-filtering.md`) or that they are unavailable
- **Output paths**: Provide paths to hotspot BED files, narrowPeak files, signal bigWig tracks, and footprint results
- **Mappability note**: Confirm which Hotspot2 mappability index was used and that it matches the read length
- **Next steps**: Suggest `motif-analysis` for TF motif enrichment in DHS peaks, or `accessibility-aggregation` for merging DHS data across samples

## For the request: "$ARGUMENTS"
