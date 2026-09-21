---
name: pipeline-hic
description: "Execute ENCODE Hi-C pipeline from FASTQ to contact matrices and loop calls. Child of pipeline-guide. Provides Nextflow execution with Docker and cloud deployment. Use when processing Hi-C data, generating contact matrices, or calling loops. Trigger on: Hi-C pipeline, chromatin conformation, contact matrix, loop calling, TAD detection, Juicer, HiCCUPS, 3D genome."
---

# ENCODE Hi-C Pipeline: FASTQ to Contact Matrices and Loops

## When to Use

- User wants to run a Hi-C processing pipeline from FASTQ to contact matrices and loop calls
- User asks about "Hi-C pipeline", "contact matrix", "loop calling", "Juicer", "HiCCUPS", or "TAD detection"
- User needs to process Hi-C data for 3D genome structure analysis
- Example queries: "process my Hi-C FASTQs", "generate contact matrices from Hi-C", "call chromatin loops with HiCCUPS"

Execute the ENCODE Hi-C pipeline for chromatin conformation capture data,
producing multi-resolution contact matrices and loop calls.

## Pipeline Overview

```
FASTQ -> FastQC (raw reads)
      -> bwa mem -SP5M (both mates in one call) -> {sample}.paired.bam
         -> pairtools parse -> sort -> dedup -> select UU
              |
              +-> Juicer pre -> .hic -> HiCCUPS -> loops (BEDPE)
              |
              +-> cooler cload + zoomify -> .mcool
```

### Not run by this workflow

Adapter trimming, TAD calling, A/B compartment calling, and every cooltools
analysis are outside this workflow. They are documented as manual, optional
steps only: trimming in `references/01-qc-trimming.md`, distance decay and
compartments in `references/04-matrix-generation.md`. cooltools and bedtools
are not installed in the container image, so those manual commands need the
conda environment (`environments/hic-env.yml` in the `bioinformatics-installer`
skill) or a separate install.

### ENCODE Repository

- **GitHub**: `ENCODE-DCC/hic-pipeline`
- **Container**: built from `scripts/Dockerfile` in this skill (`docker build -t encode-toolkit/pipeline-hic:1.0.0 scripts/`); override with `--container`
- **WDL**: Available for Cromwell execution
- **This skill**: Nextflow DSL2 reimplementation for portability

## Core Tools and Versions

Versions are those installed by `scripts/Dockerfile`, which is what the
workflow runs.

| Tool | Version | Purpose | Citation |
|------|---------|---------|----------|
| BWA-MEM | 0.7.18 | Alignment (both mates, `-SP5M`) | Li & Durbin 2009 |
| pairtools | 1.1.2 | Pair classification, dedup | Open2C |
| Juicer tools | 2.20.00 | .hic generation, HiCCUPS | Durand et al. 2016 |
| cooler | 0.9.3 | .cool/.mcool generation | Abdennur & Mirny 2020 |
| samtools | 1.19 | BAM operations | Li et al. 2009 |
| FastQC | 0.12.1 | Read quality | Andrews (Babraham) |
| MultiQC | 1.21 | Aggregated QC | Ewels et al. 2016 |

The conda alternative (`hic-env.yml`) ships no juicer_tools jar, only a JRE, so
`.hic` generation and HiCCUPS are unavailable on that route; conversely it
provides cooltools and bedtools, which the container image does not.

## Key Literature

1. **Rao et al. 2014** - "A 3D Map of the Human Genome at Kilobase Resolution
   Reveals Principles of Chromatin Looping" (Cell, ~5,000 citations)
   DOI: 10.1016/j.cell.2014.11.021

2. **Lieberman-Aiden et al. 2009** - "Comprehensive Mapping of Long-Range
   Interactions Reveals Folding Principles of the Human Genome" (Science, ~6,000 citations)
   DOI: 10.1126/science.1181369

3. **Durand et al. 2016** - "Juicer Provides a One-Click System for Analyzing
   Loop-Resolution Hi-C Experiments" (Cell Systems, ~2,000 citations)
   DOI: 10.1016/j.cels.2016.07.002

4. **Abdennur & Mirny 2020** - "Cooler: scalable storage for Hi-C data and
   other genomically labeled arrays" (Bioinformatics)
   DOI: 10.1093/bioinformatics/btz540

5. **Amemiya et al. 2019** - "The ENCODE Blacklist" (Scientific Reports, ~1,372 citations)
   DOI: 10.1038/s41598-019-45839-z

## Execution

### Quick Start (Local)

```bash
nextflow run scripts/main.nf \
    -profile local \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index '/ref/bwa_index/GRCh38.fa' \
    --chrom_sizes '/ref/hg38.chrom.sizes' \
    --outdir results/ \
    -resume
```

### SLURM HPC

```bash
nextflow run scripts/main.nf \
    -profile slurm \
    --container /path/to/pipeline-hic.sif \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index '/ref/bwa_index/GRCh38.fa' \
    --chrom_sizes '/ref/hg38.chrom.sizes' \
    --outdir results/ \
    -resume
```

### Cloud (GCP / AWS)

```bash
# Google Cloud Batch
nextflow run scripts/main.nf -profile gcp \
    --container us-docker.pkg.dev/<project>/<repo>/pipeline-hic:1.0.0 \
    --gcp_project <project> \
    --gcp_workdir gs://<bucket>/work \
    --reads 'gs://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index gs://<bucket>/ref/GRCh38.fa \
    --chrom_sizes gs://<bucket>/ref/hg38.chrom.sizes \
    --outdir gs://<bucket>/results

# AWS Batch
nextflow run scripts/main.nf -profile aws \
    --container <account>.dkr.ecr.<region>.amazonaws.com/pipeline-hic:1.0.0 \
    --aws_queue <job-queue> \
    --aws_workdir s3://<bucket>/work \
    --reads 's3://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --bwa_index s3://<bucket>/ref/GRCh38.fa \
    --chrom_sizes s3://<bucket>/ref/hg38.chrom.sizes \
    --outdir s3://<bucket>/results
```

`--outdir` only sets where results are published; Google Batch and AWS Batch
stage every task through the work directory, and the workflow stops with an
error if it or the project/queue is missing.

## Resource Requirements

| Step | CPUs | RAM | Time (2B contacts) |
|------|------|-----|---------------------|
| BWA alignment | 8 | 16 GB | 4-6 hours |
| pairtools parse + sort | 4 | 16 GB | 2-3 hours |
| pairtools dedup | 4 | 16 GB | 1-2 hours |
| Juicer pre + hic | 4 | 64 GB | 2-4 hours |
| HiCCUPS | 4 | 16 GB | 1-2 hours |
| **Total** | **8** | **64 GB** | **8-16 hours** |

The RAM column is each step's first-attempt request. Every process asks for that
much memory per attempt, so a task killed for exceeding it is retried with more
(at most two retries, capped by `--max_memory`). Failures with any other exit
status stop the run.

## Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | required | Glob pattern to paired FASTQ files |
| `--bwa_index` | required | BWA index prefix: the genome FASTA path whose `.amb .ann .bwt .pac .sa` files sit beside it (every file starting with this prefix is staged) |
| `--chrom_sizes` | required | Chromosome sizes file |
| `--outdir` | `./results` | Output directory |
| `--resolutions` | `1000,5000,10000,25000,50000,100000,250000,500000,1000000` | Matrix resolutions for `juicer_tools pre` and `cooler zoomify`: a comma-separated list of positive integers. The smallest value is the cooler base bin, and every other value must be a multiple of it; the workflow stops with an error otherwise |
| `--hiccups_resolutions` | `5000,10000,25000` | Resolutions HiCCUPS calls loops at. Only 5000, 10000 and 25000 are accepted, and each must also be listed in `--resolutions`; the workflow stops with an error otherwise. Peak width (`-p`), window width (`-i`), merge radius (`-d`) and FDR (`-f`) follow Juicer's published per-resolution defaults, one value per resolution |
| `--min_mapq` | `30` | Minimum MAPQ passed to `pairtools parse` |
| `--hiccups_gpu` | `false` | Run HiCCUPS on an NVIDIA GPU. By default the CPU mode is used, which only searches within 8 Mb of the diagonal |
| `--assembly` | `hg38` | Assembly name recorded in the `.mcool` metadata (`cooler cload --assembly`). The `.hic` file is built from the chrom.sizes file only |

### Infrastructure parameters (`nextflow.config`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--container` | `encode-toolkit/pipeline-hic:1.0.0` | Image built from `scripts/Dockerfile`. Pass a registry image for `gcp`/`aws`, or a `.sif` file for `slurm` |
| `--max_cpus`, `--max_memory`, `--max_time` | `16`, `128.GB`, `48.h` | Upper bounds applied to every process |
| `--slurm_queue`, `--slurm_account` | `normal`, none | SLURM partition and account |
| `--gcp_project`, `--gcp_workdir` | none (both required for `-profile gcp`) | Google Cloud project and `gs://` work directory |
| `--gcp_location`, `--gcp_disk` | `us-central1`, `200.GB` | Google Batch region and per-task disk |
| `--aws_queue`, `--aws_workdir` | none (both required for `-profile aws`) | AWS Batch job queue and `s3://` work directory |
| `--aws_region`, `--aws_cli_path` | `us-east-1`, `/home/ec2-user/miniconda/bin/aws` | AWS region, and the AWS CLI path inside the Batch AMI |

## Output Files

```
results/
  fastqc/
    *_fastqc.html                 # Raw read quality
    *_fastqc.zip
  alignment/
    {sample}.paired.bam           # Both mates from one bwa mem -SP5M call
  pairs/
    {sample}.parse_stats.txt      # pairtools parse stats (pair-type breakdown)
    {sample}.dedup.pairs.gz       # Classified, sorted, deduplicated pairs
    {sample}.dedup_stats.txt      # pairtools dedup stats (duplication, complexity)
  matrices/
    {sample}.hic                  # Juicer .hic file (primary output)
    {sample}.mcool                # Cooler multi-resolution matrix
  loops/
    {sample}.hiccups_loops.bedpe  # HiCCUPS merged_loops.bedpe, renamed
  qc/
    {sample}.contact_stats.txt    # pairtools stats on the selected UU pairs
  multiqc/
    multiqc_report.html
  pipeline_info/
    timeline.html
    report.html
    trace.txt
```

The UU-selected pairs file and the `pairtools stats` run on it are intermediate:
only `qc/{sample}.contact_stats.txt` is published, not the selected pairs.

### .hic File Format

The .hic format (Juicer) stores multi-resolution contact matrices with
normalization vectors. Can be visualized in Juicebox and loaded by
`hic-straw` in Python/R.

### .mcool File Format

The .mcool format (cooler) is an HDF5-based multi-resolution contact matrix.
Widely supported by `cooler`, `cooltools`, `HiGlass`, and `FAN-C`.

## QC Thresholds (ENCODE Standards)

This is the only QC threshold table for this skill; the reference files point
back to it.

| Metric | Pass | Warning | Fail | Computed from |
|--------|------|---------|------|---------------|
| Valid (UU) pair fraction | >40% | 25-40% | <25% | `pairs/{sample}.parse_stats.txt` |
| Cis contacts (>20kb) | >40% | 25-40% | <25% | `qc/{sample}.contact_stats.txt` |
| Cis/trans ratio | >1.5 | 1.0-1.5 | <1.0 | `qc/{sample}.contact_stats.txt` |
| Library complexity (unique/total) | >0.7 | 0.5-0.7 | <0.5 | `pairs/{sample}.dedup_stats.txt` |

`qc/{sample}.contact_stats.txt` is computed after UU selection, so its pair-type
breakdown is 100% UU by construction. Read pair types from
`pairs/{sample}.parse_stats.txt` instead.

### Resolution vs Depth Requirements

| Resolution | Minimum Contacts Needed | Typical Depth |
|------------|------------------------|---------------|
| 1 kb | >2 billion | Very deep |
| 5 kb | >500 million | Deep |
| 10 kb | >200 million | Standard |
| 25 kb | >50 million | Moderate |
| 100 kb | >10 million | Low |

## Pair Classification

pairtools assigns each read pair a two-letter code (one letter per side:
U unique, R rescued, M multi, N null/unmapped, W walk, D duplicate, X corrupt):

| Category | Description | Use |
|----------|-------------|-----|
| UU | Both sides uniquely mapped | Valid contact -- the only type this workflow keeps |
| UR / RU | One unique, one rescued | Valid but not selected here |
| NU | One unique, one unmapped | Not used |
| NM | One unmapped, one multi-mapped | Not used |
| MM | Both multi-mapped | Not used |
| WW | Complex walk (multiple ligation events), masked by `--walks-policy mask` | Not used |
| DD | Duplicate | Removed by `pairtools dedup` |
| XX | Corrupt record | Not used |

This workflow selects UU only
(`pairtools select '(pair_type == "UU")'`) before matrix generation.

## Critical Pitfalls

### Restriction Enzyme Choice
The restriction enzyme determines fragment size and resolution:
- **MboI/DpnII** (GATC): 4-cutter, ~256 bp average fragment -- higher resolution
- **HindIII** (AAGCTT): 6-cutter, ~4 kb average fragment -- lower resolution
- **Arima** (proprietary): Two enzymes, ~160 bp average -- highest resolution
- Always verify which enzyme was used before interpreting resolution
- The workflow itself is enzyme-agnostic: pairtools works at read-pair level and the `.hic`
  file is built without a restriction-site file, so there is no enzyme parameter to set

### Normalization Method
Different normalization methods yield different results:
- **KR** (Knight-Ruiz): built by `juicer_tools pre -k KR,VC,VC_SQRT` and used by HiCCUPS (`-k KR`)
- **ICE** (Imakaev et al.): applied to the .mcool by `cooler zoomify --balance`
- **VC** (Vanilla Coverage): simple coverage normalization, also built into the .hic
- Always document which normalization a downstream analysis read.

### Resolution Depends on Depth
Do not call features at resolutions unsupported by sequencing depth:
- Calling 1 kb loops from 100M contacts will produce noise
- Check the Juicer resolution QC to determine achievable resolution
- HiCCUPS runs at the resolutions in `--hiccups_resolutions` (default 5 kb, 10 kb
  and 25 kb); each of them must also be listed in `--resolutions`, because HiCCUPS
  reads them out of the `.hic` file

### Ligation Artifacts
Monitor the pair-type breakdown in `pairs/{sample}.parse_stats.txt`:
- WW pairs are complex walks (more than one ligation in a read); `--walks-policy mask`
  masks them so they never reach the contact matrix
- A large unmapped/multi-mapped fraction points at poor library or the wrong genome
- The workflow produces no re-ligation distance plot; derive one manually from the
  published pairs file if needed

## Provenance Integration

After pipeline completion, log all outputs:

```python
encode_log_derived_file(
    file_path="/results/matrices/sample1.hic",
    source_accessions=["ENCSR...", "ENCFF..."],
    description="Hi-C contact matrix from ENCODE Hi-C pipeline",
    file_type="hic",
    tool_used="BWA 0.7.18 + pairtools 1.1.2 + Juicer 2.20.00",
    parameters="--min_mapq 30, UU pairs only, KR/VC/VC_SQRT normalization, resolutions 1kb-1Mb"
)
```

## Reference Files

Detailed step-by-step documentation is provided in the `references/` directory:

1. `01-qc-trimming.md` -- Read QC (trimming is a manual option, not run here)
2. `02-alignment.md` -- BWA `-SP5M` alignment of both mates in one call
3. `03-pair-processing.md` -- pairtools parse, sort, dedup, and select
4. `04-matrix-generation.md` -- Juicer .hic and cooler .mcool generation; manual cooltools analyses
5. `05-loop-calling.md` -- HiCCUPS loop detection and QC

## Walkthrough: Processing ENCODE Hi-C from FASTQ to Contact Maps and Loops

**Goal**: Process raw Hi-C FASTQ files through the ENCODE pipeline to generate contact matrices and chromatin loop calls.
**Context**: Hi-C captures 3D chromatin organization. The pipeline uses BWA for chimeric read alignment, pairtools for pair processing, and Juicer/HiCCUPS for loop calling.

### Step 1: Find Hi-C experiment

```
encode_get_experiment(accession="ENCSR000AKA")
```

Expected output:
```json
{
  "accession": "ENCSR000AKA",
  "assay_title": "Hi-C",
  "biosample_summary": "GM12878",
  "bio_replicate_count": 2,
  "status": "released"
}
```

### Step 2: List FASTQ files

```
encode_list_files(experiment_accession="ENCSR000AKA", file_format="fastq")
```

Expected output (a JSON array of file records; fields abridged):
```json
[
  {"accession": "ENCFF500HI1", "file_format": "fastq", "output_type": "reads", "biological_replicates": [1], "file_size_human": "34.2 GB", "status": "released"},
  {"accession": "ENCFF501HI2", "file_format": "fastq", "output_type": "reads", "biological_replicates": [1], "file_size_human": "35.2 GB", "status": "released"}
]
```

**Interpretation**: Hi-C paired-end reads represent chimeric ligation junctions. Each read pair captures a 3D contact.

### Step 3: Name the files so a read-pair glob can find them

ENCODE FASTQs are named by accession, so the two mates of a pair share no
prefix, and the workflow matches file pairs with a `{1,2}` glob. Which mate a
file is comes from its page on encodeproject.org (`paired_end` 1 or 2, and
`paired_with` naming the other accession), not from any tool here. Link the
files into the shape the glob expects:

```bash
mkdir -p fastq
ln -s "$PWD/ENCFF500HI1.fastq.gz" fastq/ENCSR000AKA_R1.fastq.gz
ln -s "$PWD/ENCFF501HI2.fastq.gz" fastq/ENCSR000AKA_R2.fastq.gz
```

### Step 4: Run the Hi-C pipeline

```bash
nextflow run scripts/main.nf \
  -profile local \
  --reads 'fastq/ENCSR000AKA_R{1,2}.fastq.gz' \
  --bwa_index '/ref/bwa_index/GRCh38.fa' \
  --chrom_sizes '/ref/hg38.chrom.sizes' \
  --outdir results/ \
  -resume
```

Key pipeline steps:
1. FastQC on the raw reads
2. BWA-MEM `-SP5M` alignment of both mates in one call (chimeric read handling)
3. pairtools parse + sort (classify pairs, MAPQ 30, mask walks)
4. pairtools dedup (remove PCR duplicates), then select UU pairs
5. Contact matrix generation (`.hic` via Juicer, `.mcool` via cooler)
6. Loop calling (HiCCUPS at the `--hiccups_resolutions`, by default 5 kb, 10 kb
   and 25 kb, merged into one BEDPE)

### Step 5: Validate output quality

Use the QC threshold table above with `pairs/{sample}.parse_stats.txt`,
`pairs/{sample}.dedup_stats.txt` and `qc/{sample}.contact_stats.txt`.

### Step 6: Identify significant loops

Download loop calls for downstream analysis:
```
encode_list_files(experiment_accession="ENCSR000AKA", file_format="bedpe", assembly="GRCh38")
```

### Integration with downstream skills
- Loop calls (BEDPE) feed into -> **hic-aggregation** for cross-tissue loop catalog
- Loop anchors feed into -> **peak-annotation** for enhancer-promoter assignment
- Contact data integrates with -> **visualization-workflow** for 3D genome display
- Pipeline provenance logged by -> **data-provenance**

## Code Examples

### 1. Find Hi-C data for 3D genome analysis

```
encode_search_experiments(
  assay_title="Hi-C",
  organ="heart"
)
```

Expected output:
```json
{
  "results": [
    {
      "accession": "ENCSR654HRT",
      "assay_title": "Hi-C",
      "biosample_summary": "heart left ventricle tissue male adult (51 years)",
      "status": "released"
    }
  ],
  "total": 4,
  "limit": 25,
  "offset": 0,
  "has_more": false,
  "next_offset": null
}
```

### 2. Get experiment details for pipeline configuration

```
encode_get_experiment(accession="ENCSR654HRT")
```

Expected output:
```json
{
  "accession": "ENCSR654HRT",
  "assay_title": "Hi-C",
  "bio_replicate_count": 2,
  "biosample_summary": "heart left ventricle tissue male adult (51 years)",
  "assembly": ["GRCh38"],
  "audit_error_count": 0,
  "audit_warning_count": 1
}
```

## Integration

| This skill produces... | Feed into... | Purpose |
|---|---|---|
| Chromatin loops (BEDPE) | **hic-aggregation** | Cross-tissue loop catalog |
| Loop anchors (BED) | **peak-annotation** | Assign genes to loop-connected enhancers |
| Contact matrices (.hic / .mcool) | **visualization-workflow** | 3D genome visualization |
| Loop-disrupting coordinates | **variant-annotation** | Identify variants breaking chromatin contacts |
| QC metrics | **quality-assessment** | Validate Hi-C library quality |
| Pipeline parameters | **data-provenance** | Record BWA/pairtools/Juicer versions |
| Loop anchor regions | **motif-analysis** | Discover CTCF motifs at loop anchors |

## Related Skills

- `pipeline-guide` -- Parent skill with compute resource assessment and cloud setup
- `hic-aggregation` -- Aggregate Hi-C loops across samples/tissues
- `quality-assessment` -- Evaluate pipeline output quality metrics
- `data-provenance` -- Track all pipeline inputs, outputs, and parameters
- `download-encode` -- Download ENCODE Hi-C FASTQ files for pipeline input
- `publication-trust` -- Verify literature claims backing analytical decisions

## Presenting Results

When reporting Hi-C pipeline results:

- **Valid pair count**: Report the UU pair count and its fraction of all parsed pairs from `pairs/{sample}.parse_stats.txt`. UU is the only pair type this workflow carries forward
- **Cis/trans ratio**: Report the cis/trans contact ratio (>1.5 pass) and long-range cis fraction (>20kb, >40% expected) from `qc/{sample}.contact_stats.txt`. These are the primary Hi-C quality indicators
- **Contact matrix resolution**: Report the achievable resolution based on sequencing depth (e.g., "500M valid pairs supports 5kb resolution") and list the `--resolutions` actually generated
- **Loop counts**: Report the number of loops in `loops/{sample}.hiccups_loops.bedpe` (the file starts with a `#chr1 ...` header line, so exclude it from the count). HiCCUPS searches the resolutions in `--hiccups_resolutions` (default 5 kb, 10 kb and 25 kb) and merges them into that single file; per-resolution files stay in the Nextflow work directory
- **Matrix paths**: Provide paths to the .hic file (Juicebox-compatible) and .mcool file (cooler/HiGlass-compatible)
- **Key QC metrics**: Present library complexity (unique/total >0.7, from `pairs/{sample}.dedup_stats.txt`) and the pair-type breakdown from `pairs/{sample}.parse_stats.txt` in a summary table
- **Normalization**: Note that `.hic` carries KR, VC and VC_SQRT vectors (HiCCUPS uses KR) and the `.mcool` is ICE-balanced by `cooler zoomify --balance`
- **Not produced here**: TAD calls, A/B compartments and cooltools outputs are not generated by this workflow; say so rather than implying they are missing
- **Next steps**: Suggest `hic-aggregation` for cross-sample loop catalogs, or `visualization-workflow` for Juicebox/HiGlass session setup

## For the request: "$ARGUMENTS"
