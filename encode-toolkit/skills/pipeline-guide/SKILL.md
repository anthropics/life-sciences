---
name: pipeline-guide
description: Access ENCODE uniform analysis pipelines, generate user-specific Nextflow/WDL pipelines, manage compute resources, and integrate with cloud platforms. Use when the user wants to understand ENCODE pipelines, run pipelines on their own data, generate custom Nextflow workflows from ENCODE pipeline code, check compute requirements (CPU/GPU/memory), run pipelines in background, or integrate with Google Cloud, AWS, or other cloud platforms. Also use when the user asks about ENCODE pipeline outputs, processing standards, software versions, or wants to replicate ENCODE processing. Covers local execution, HPC, and cloud deployment with resource-aware scheduling. Use this skill for ANY pipeline execution, workflow generation, or compute resource management task involving ENCODE data.
---

# ENCODE Pipeline Guide and Custom Workflow Generation

## When to Use

- User wants to understand ENCODE uniform analysis pipelines or run them on their own data
- User asks about "ENCODE pipeline", "Nextflow", "WDL", "processing standards", or "pipeline requirements"
- User needs to generate a custom Nextflow/WDL workflow based on ENCODE pipeline specifications
- User wants to know compute requirements (CPU, GPU, memory, storage) for running pipelines
- Example queries: "how do I run the ENCODE ChIP-seq pipeline?", "what are the compute requirements for Hi-C processing?", "generate a Nextflow pipeline for my ATAC-seq data"

Understand ENCODE pipelines, generate user-specific workflows in Nextflow/WDL, and manage compute resources for local, HPC, and cloud execution.

## ENCODE Uniform Analysis Pipelines

ENCODE uses standardized pipelines for each assay type, ensuring reproducibility across all datasets. All pipelines are:
- **Open source**: GitHub (github.com/ENCODE-DCC)
- **Containerized**: Docker and Singularity images
- **Written in WDL**: Workflow Description Language (Cromwell execution engine)
- **Portable**: Local, HPC (SLURM, SGE, PBS), or cloud (Google Cloud, AWS, Azure)

### Pipeline Repository Map

The official ENCODE pipelines are the WDL workflows below. The `pipeline-*` skills in this
toolkit are independent Nextflow implementations that follow the same standards; they are
written and maintained by the ENCODE Toolkit author, not by the ENCODE DCC. Each skill builds
its own image from the `scripts/Dockerfile` it ships with.

| Assay | Official ENCODE pipeline (WDL) | Primary Tools | Toolkit skill (Nextflow) |
|-------|-------------------------------|---------------|--------------------------|
| ChIP-seq | `ENCODE-DCC/chip-seq-pipeline2` | BWA, MACS2, IDR | `pipeline-chipseq` |
| ATAC-seq | `ENCODE-DCC/atac-seq-pipeline` | Bowtie2, MACS2, IDR | `pipeline-atacseq` |
| RNA-seq | `ENCODE-DCC/rna-seq-pipeline` | STAR, RSEM | `pipeline-rnaseq` |
| DNase-seq | `ENCODE-DCC/dnase-seq-pipeline` | BWA, Hotspot2 | `pipeline-dnaseseq` |
| WGBS | `ENCODE-DCC/dna-me-pipeline` | Bismark, MethylDackel | `pipeline-wgbs` |
| Hi-C | `ENCODE-DCC/hic-pipeline` | BWA, Juicer, HiCCUPS | `pipeline-hic` |
| CUT&RUN | none published by ENCODE | Bowtie2, SEACR/MACS2 | `pipeline-cutandrun` |
| scRNA-seq | see the ENCODE portal pipeline pages | STARsolo | — |
| scATAC-seq | see the ENCODE portal pipeline pages | Chromap | — |

ENCODE also publishes images for some pipelines on Docker Hub, for example
`encodedcc/chip-seq-pipeline:v2.2.1` and `encodedcc/atac-seq-pipeline:v2.2.0`. Those images are
built for the WDL workflows and are not what the Nextflow skills here run.

## Literature Foundation

| Reference | Year | Relevance | Citations |
|-----------|------|-----------|-----------|
| Di Tommaso et al. "Nextflow enables reproducible computational workflows" | 2017 | Nextflow workflow manager | ~2,800 |
| Ewels et al. "The nf-core framework for community-curated bioinformatics pipelines" | 2020 | nf-core community pipelines | ~1,900 |
| Kurtzer et al. "Singularity: Scientific containers for mobility of compute" | 2017 | Singularity containers for HPC | ~2,500 |
| Merkel "Docker: lightweight Linux containers for consistent development and deployment" | 2014 | Docker containerization | ~3,000 |
| ENCODE Project Consortium "Expanded encyclopaedias of DNA elements" | 2020 | ENCODE Phase 3 standards | ~1,200 |
| Gruening et al. "Bioconda: sustainable and comprehensive software distribution" | 2018 | Bioconda packaging ecosystem | ~1,400 |

## Pipeline Output Types by Assay

### ChIP-seq Pipeline
| Output Type | Format | Description | Use For |
|------------|--------|-------------|---------|
| alignments | bam | Filtered, deduplicated | Reprocessing, visualization |
| signal of unique reads | bigWig | Unique read signal | Genome browser |
| fold change over control | bigWig | Normalized signal | Comparative visualization |
| IDR thresholded peaks | bed narrowPeak | Reproducible peaks | Peak analysis (gold standard) |
| pseudoreplicated peaks | bed narrowPeak | Single-replicate peaks | When only 1 replicate |
| optimal IDR peaks | bed narrowPeak | Pooled replicate peaks | Most complete peak set |

### ATAC-seq Pipeline
| Output Type | Format | Description | Use For |
|------------|--------|-------------|---------|
| alignments | bam | No-mito, deduplicated | Reprocessing |
| signal of unique reads | bigWig | Signal track | Genome browser |
| IDR thresholded peaks | bed narrowPeak | Reproducible peaks | Accessibility analysis |
| pseudoreplicated peaks | bed narrowPeak | Single-replicate | Backup peaks |

### RNA-seq Pipeline
| Output Type | Format | Description | Use For |
|------------|--------|-------------|---------|
| alignments | bam | STAR-aligned | Visualization, reprocessing |
| gene quantifications | tsv | Gene-level counts (RSEM) | Differential expression |
| transcript quantifications | tsv | Transcript-level counts | Isoform analysis |
| signal of unique reads | bigWig | Strand-specific signal | Genome browser |

### WGBS Pipeline
| Output Type | Format | Description | Use For |
|------------|--------|-------------|---------|
| alignments | bam | Bisulfite-converted | Reprocessing |
| methylation state at CpG | bed bedMethyl | Per-CpG levels | Methylation analysis |

### Hi-C Pipeline
| Output Type | Format | Description | Use For |
|------------|--------|-------------|---------|
| contact matrix | hic | Interaction frequencies | TAD/compartment calling |
| loops | bedpe | Called loops | Loop analysis |

## Choosing the Right Output Files

### Decision Table
| Analysis Goal | File Type | Output Type | Priority |
|--------------|-----------|-------------|----------|
| Visualization | bigWig | fold change over control (ChIP) / signal of unique reads (others) | preferred_default=True |
| Peak overlap | bed narrowPeak | IDR thresholded peaks | Highest confidence |
| Quantitative | tsv / bed | gene quantifications / methylation state | Pipeline defaults |
| Custom processing | fastq | reads | When ENCODE pipeline doesn't match |

```
encode_list_files(experiment_accession="ENCSR...", preferred_default=True)
```

## Step 1: Assess User Compute Resources

Before generating any pipeline, check available resources:

### System Check Commands
```bash
# CPU cores
nproc                              # Linux
sysctl -n hw.ncpu                  # macOS

# Memory
free -h                            # Linux
sysctl -n hw.memsize | awk '{print $1/1024/1024/1024 " GB"}'  # macOS

# Disk space
df -h /path/to/data/

# GPU (if applicable)
nvidia-smi                         # NVIDIA GPU
# Note: Most ENCODE pipelines do NOT require GPU

# Docker availability
docker --version
docker info | grep "Total Memory"

# Singularity (for HPC)
singularity --version
```

### Minimum Resource Requirements by Pipeline

| Pipeline | Min CPU | Min RAM | Min Disk | GPU | Time Estimate (per sample) |
|----------|---------|---------|----------|-----|---------------------------|
| ChIP-seq | 4 cores | 16 GB | 50 GB | No | 2–4 hours |
| ATAC-seq | 4 cores | 16 GB | 50 GB | No | 2–4 hours |
| RNA-seq | 8 cores | 32 GB | 100 GB | No | 4–8 hours (index build) |
| WGBS | 8 cores | 48 GB | 200 GB | No | 12–24 hours |
| Hi-C | 8 cores | 64 GB | 200 GB | No | 8–16 hours |
| DNase-seq | 8 cores | 16 GB | see `pipeline-dnaseseq` | No | 3–6 hours |
| CUT&RUN | 8 cores | 8 GB | see `pipeline-cutandrun` | No | 1.5–3 hours |

The DNase-seq and CUT&RUN figures are the per-sample totals from those skills' own resource
tables. There is no toolkit pipeline skill for scRNA-seq or scATAC-seq, so no row is given.

### Resource Scaling
- **CPU**: Alignment steps are parallelizable; doubling cores approximately halves alignment time
- **RAM**: Genome index loading is the bottleneck; STAR requires ~32 GB for human genome
- **Disk**: FASTQ + BAM + intermediate files can exceed 100 GB per sample
- **Network**: ENCODE downloads at ~50–200 MB/s; plan for transfer time

## Step 2: Generate Custom Nextflow Workflows

When the user needs to run ENCODE-style processing, generate Nextflow workflows that mirror ENCODE pipeline logic.

### Why Nextflow Over WDL
- **Broader adoption**: Nextflow is used by nf-core, most HPC centers, and cloud platforms
- **Native container support**: Docker, Singularity, Podman
- **Cloud integration**: AWS Batch, Google Cloud Batch, Azure Batch natively
- **Resource management**: Built-in CPU/memory/time limits per process
- **Resume capability**: Failed runs restart from last successful step

### Nextflow Pipeline Template

```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Pipeline parameters
params.reads         = null          // Input FASTQ glob, e.g. '*_R{1,2}.fastq.gz'
params.genome        = 'GRCh38'      // Genome assembly
params.bwa_index     = null          // BWA index prefix or directory
params.outdir        = './results'   // Output directory

// Example: ChIP-seq alignment process
process ALIGN_READS {
    tag "${sample_id}"
    cpus 8
    memory { 16.GB * task.attempt }
    time   { 4.h * task.attempt }
    container params.container

    input:
    tuple val(sample_id), path(reads)
    path genome_index

    output:
    tuple val(sample_id), path("*.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} -M ${genome_index}/${params.genome}.fa ${reads} | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}
```

Per-process resources go on the process (or in `withName:` blocks in the config); the global
ceiling belongs in `nextflow.config` via `process.resourceLimits`, shown next.

### Resource-Aware Configuration

Generate a `nextflow.config` based on the user's system. This mirrors the `nextflow.config` that
every `pipeline-*` skill ships; the profile names are `local`, `slurm`, `gcp` and `aws` — no
others exist. `params.container` used by the process above is declared here.

```nextflow
params {
    outdir     = './results'

    // Auto-detected from the user's system
    max_cpus   = ${detected_cpus}
    max_memory = '${detected_memory}.GB'
    max_time   = '72.h'

    // Image built from the skill's scripts/Dockerfile. Override with a registry image for
    // gcp/aws, or a .sif file for slurm.
    container  = 'encode-toolkit/pipeline-chipseq:1.0.0'

    // Scheduler and cloud settings (only read by the matching profile)
    slurm_queue   = 'normal'
    slurm_account = null
    gcp_project   = null
    gcp_location  = 'us-central1'
    gcp_workdir   = null
    gcp_disk      = '200.GB'
    aws_queue     = null
    aws_region    = 'us-east-1'
    aws_workdir   = null
    aws_cli_path  = '/home/ec2-user/miniconda/bin/aws'
}

process {
    container      = params.container
    // Retry only when a task was killed for exceeding its memory or time (exit codes
    // 130-145 and 104); every process requests memory per attempt, so the retry gets more.
    // Any other failure is a real error: stop and report it.
    errorStrategy  = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries     = 2
    // Nextflow's built-in ceiling: every process request is capped at these values.
    resourceLimits = [cpus: params.max_cpus, memory: params.max_memory, time: params.max_time]
}

profiles {
    local {
        process.executor  = 'local'
        docker.enabled    = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    slurm {
        process.executor       = 'slurm'
        process.queue          = params.slurm_queue
        process.clusterOptions = params.slurm_account ? "--account=${params.slurm_account}" : null
        singularity.enabled    = true
        singularity.autoMounts = true
    }

    gcp {
        process.executor  = 'google-batch'
        process.disk      = params.gcp_disk
        google.project    = params.gcp_project
        google.location   = params.gcp_location
        google.batch.spot = true
        workDir           = params.gcp_workdir
    }

    aws {
        process.executor  = 'awsbatch'
        process.queue     = params.aws_queue
        aws.region        = params.aws_region
        aws.batch.cliPath = params.aws_cli_path
        workDir           = params.aws_workdir
    }
}
```

`process.resourceLimits` replaces the hand-written `check_max()` helper that older nf-core
configs use: it applies to `cpus`, `memory` and `time` alike, so there is no helper to keep in
sync.

Google Batch and AWS Batch stage every task through object storage, so `-profile gcp` needs
`--gcp_project` and `--gcp_workdir gs://<bucket>/work`, and `-profile aws` needs `--aws_queue`
and `--aws_workdir s3://<bucket>/work`. `--outdir` only sets where results are published. Cloud
runs also need `--container <registry image>`; the default image name is local-only.

## Step 3: Cloud Integration

### Available Integrations (Official Marketplace)

For users who cannot run pipelines locally, offer cloud integration:

#### Google Cloud / Colab
- **Nextflow + Google Cloud Batch**: Run full pipelines on Google Cloud (`-profile gcp`)
- **Google Colab**: For interactive analysis (R/Python notebooks)
  - Limited to 12 GB RAM (free tier) or 25 GB (Pro)
  - GPU available (useful for deep learning, not standard pipelines)
  - Best for: downstream analysis after pipeline completion

#### AWS
- **Nextflow + AWS Batch**: Run pipelines on AWS (`-profile aws`)
- **AWS SageMaker**: For ML-based analysis
- Best for: Large-scale batch processing

#### Other Platforms
- **Terra (Broad Institute)**: WDL-native platform, ENCODE pipelines pre-installed
- **DNAnexus**: Cloud genomics platform with ENCODE pipeline apps
- **Galaxy**: Web-based, no coding required

### Cloud Cost Estimates
| Pipeline | Cloud Instance | Estimated Cost/Sample |
|----------|---------------|---------------------|
| ChIP-seq | n1-standard-8 (GCP) / m5.2xlarge (AWS) | $2–5 |
| ATAC-seq | n1-standard-8 / m5.2xlarge | $2–5 |
| RNA-seq | n1-standard-16 / m5.4xlarge | $5–10 |
| WGBS | n1-highmem-16 / r5.4xlarge | $10–25 |
| Hi-C | n1-highmem-16 / r5.4xlarge | $8–20 |

DNase-seq and CUT&RUN are not listed: no cost measurements exist for them. Estimate from their
resource rows above (both fit an 8-core instance) and your provider's current rates; the `gcp`
profile enables Batch spot instances by default.

## Step 4: Background Execution

### Local Background Execution
```bash
# Nextflow's own -bg flag detaches the run; redirect its log and you do not need nohup.
# Pass every required parameter for the pipeline you are running (see its SKILL.md).
nextflow run pipeline-chipseq/scripts/main.nf \
    -profile local \
    --reads '/path/to/reads/*_R{1,2}.fastq.gz' \
    --chrom_sizes /ref/hg38.chrom.sizes \
    --outdir results/ \
    -resume \
    -bg \
    > pipeline.log 2>&1

# Monitor progress
tail -f pipeline.log
nextflow log last
```

### Screen/tmux for Long Runs
```bash
# Create a persistent session
screen -S encode_pipeline
# or
tmux new -s encode_pipeline

# Run pipeline inside session (pass every required parameter for that pipeline)
nextflow run pipeline-chipseq/scripts/main.nf -profile local \
    --reads '/path/to/reads/*_R{1,2}.fastq.gz' \
    --chrom_sizes /ref/hg38.chrom.sizes \
    --outdir results/ -resume

# Detach: Ctrl+A then D (screen) or Ctrl+B then D (tmux)
# Reattach later: screen -r encode_pipeline / tmux attach -t encode_pipeline
```

## Step 5: Extract ENCODE Pipeline Code Snippets

When the user needs specific processing steps (not full pipelines), extract the relevant code:

### Common Snippets

#### Alignment (ChIP-seq / ATAC-seq)
```bash
# What pipeline-chipseq runs (BWA_MEM). The -F 1804 flag filter is a separate
# step (FILTER_SORT), not part of this command.
bwa mem -t ${NCPUS} ${GENOME_INDEX}/${GENOME}.fa ${FASTQ_R1} ${FASTQ_R2} | \
    samtools view -@ ${NCPUS} -bS -q 30 - | \
    samtools sort -@ ${NCPUS} -m 2G -o sample.bam -
samtools index sample.bam
samtools flagstat sample.bam > sample.flagstat.txt

# What pipeline-atacseq runs (BOWTIE2_ALIGN)
bowtie2 --very-sensitive -X 2000 --no-mixed --no-discordant \
    --threads ${NCPUS} -x ${GENOME_INDEX}/${GENOME} \
    -1 ${FASTQ_R1} -2 ${FASTQ_R2} 2> sample.bowtie2.log | \
    samtools view -@ ${NCPUS} -bS -q 30 -f 2 - | \
    samtools sort -@ ${NCPUS} -m 2G -o sample.bam -
samtools index sample.bam

# Mark/remove duplicates (MARK_DUPLICATES in both workflows)
picard MarkDuplicates \
    INPUT=sample.filtered.bam \
    OUTPUT=sample.dedup.bam \
    METRICS_FILE=sample.dup_metrics.txt \
    REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
```

These are the commands the toolkit's own Nextflow workflows run
(`pipeline-chipseq/scripts/main.nf`, `pipeline-atacseq/scripts/main.nf`). They
are not copied from ENCODE's WDL pipelines, whose commands differ.

#### Peak Calling (MACS2)
```bash
# What pipeline-chipseq runs. -f is BAMPE for paired-end, BAM for single-end;
# -g is 'hs' (GRCh38) or 'mm' (mm10).
macs2 callpeak \
    -t treatment.bam -c control.bam \
    -f BAMPE -g hs -n sample \
    --qvalue 0.05 --nomodel --keep-dup all \
    --call-summits -B

# Broad marks (--peak_type broad) swap --call-summits for:
#   --broad --broad-cutoff 0.1
```

`--shift`/`--extsize` are the ATAC/single-end recipe and are not used here: they have no effect
in `BAMPE` mode, where MACS2 takes the fragment from the read pair. `pipeline-atacseq` also runs
`-f BAMPE` without them, because it applies the Tn5 offset upstream with
`alignmentSieve --ATACshift`.

#### IDR Analysis
```bash
# What pipeline-chipseq and pipeline-atacseq run, once for every pair of samples
# matched by --reads (narrow peaks only; 2 samples -> 1 comparison, 3 -> 3).
idr --samples rep1_peaks.narrowPeak rep2_peaks.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file rep1_vs_rep2.idr_peaks.txt \
    --plot \
    --idr-threshold 0.05
```

Each comparison is published as `peaks/idr/<sampleA>_vs_<sampleB>.idr_peaks.txt` (plus a
`.png`), with the two names in alphabetical order; with a single sample IDR is skipped. There is
no pooled or pseudoreplicate analysis and no rescue/self-consistency ratio in these workflows;
add those steps yourself if you need the full ENCODE IDR protocol.

#### RNA-seq Quantification
```bash
# What pipeline-rnaseq runs (STAR 2-pass + RSEM)
STAR --genomeDir ${STAR_INDEX} \
    --readFilesIn ${FASTQ_R1} ${FASTQ_R2} \
    --readFilesCommand zcat \
    --runThreadN ${NCPUS} \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --quantMode TranscriptomeSAM GeneCounts \
    --twopassMode Basic \
    --outWigType bedGraph \
    --outWigStrand Stranded \
    --outFileNamePrefix sample.

rsem-calculate-expression \
    --paired-end \
    --bam \
    --no-bam-output \
    --estimate-rspd \
    --strandedness reverse \
    --num-threads ${NCPUS} \
    sample.Aligned.toTranscriptome.out.bam \
    ${RSEM_INDEX} \
    sample
```

`--twopassMode Basic` is what makes this "STAR 2-pass". `--strandedness` must match the library
(`reverse` for dUTP protocols, `forward`, or `none`); `--outWigStrand` becomes `Unstranded` when
it is `none`. The annotation is baked into the STAR and RSEM indexes when they are built, so
there is no GTF argument here.

#### Liftover (GRCh37 → GRCh38)
```bash
# Download chain file
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Run liftover
liftOver input_hg19.bed hg19ToHg38.over.chain.gz output_hg38.bed unmapped.bed

# Log: liftOver version (Kent et al. 2002, Genome Research)
# Log: chain file source and date accessed
# Log: input count, output count, unmapped count
```

## Step 6: Language-Specific Integration

### R / Bioconductor
For users working in R, ENCODE data integrates with:
```r
# Key Bioconductor packages for ENCODE data
library(GenomicRanges)      # Genomic intervals
library(rtracklayer)        # Import BED/bigWig
library(DESeq2)             # Differential expression
library(DiffBind)           # Differential binding (ChIP-seq)
library(ChIPseeker)         # Peak annotation
library(chromVAR)           # Chromatin accessibility
library(BSgenome.Hsapiens.UCSC.hg38)  # Genome sequence
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Gene models

# Import ENCODE peak file
peaks <- rtracklayer::import("ENCFF123ABC.bed", format="narrowPeak")

# Import ENCODE bigWig signal
signal <- rtracklayer::import("ENCFF456DEF.bigWig", format="bigWig")
```

Check package availability:
```r
# CRAN
available.packages(repos="https://cran.r-project.org")[,"Version"]

# Bioconductor
BiocManager::available()
BiocManager::version()
```

### Python
```python
# Key Python packages for ENCODE data
import pyBigWig          # Read bigWig files
import pybedtools        # BED operations
import pysam             # BAM file access
import scanpy as sc      # Single-cell analysis
import anndata           # AnnData format
import cooler            # Hi-C contact matrices
import pydeseq2          # Differential expression

# Import ENCODE peak file
import pandas as pd
peaks = pd.read_csv("ENCFF123ABC.bed", sep="\t", header=None,
                     names=["chr","start","end","name","score","strand",
                            "signalValue","pValue","qValue","peak"])
```

### Bash / Command Line
Core tools for ENCODE data processing:
```bash
# Essential tools and typical versions
bedtools --version    # v2.31.0 - genomic arithmetic
samtools --version    # 1.19 - BAM/CRAM operations
tabix                 # indexing BED/VCF
bigWigToBedGraph      # UCSC Kent tools
bedToBigBed           # UCSC Kent tools
macs2 --version       # 2.2.9.1 - peak calling
idr --version         # 2.0.4.2 - reproducibility
deeptools --version   # 3.5.5 - signal visualization
```

## Provenance Integration

When generating or running any pipeline, integrate with the data-provenance skill:

1. **Before execution**: Log all input files, tool versions, reference files
2. **During execution**: Capture stdout/stderr, resource usage
3. **After execution**: Log all output files with MD5 checksums, record runtime
4. **Script storage**: Save the generated pipeline script in `scripts/` directory

Every pipeline run should produce a provenance entry that enables methods writing.

## Pitfalls and Edge Cases

### Version Mismatches
- ENCODE has used multiple pipeline versions over the years
- Files from different pipeline versions may not be directly comparable
- Check the `analysis` field in file metadata for pipeline version
- When reprocessing, use the same pipeline version as ENCODE for comparability

### Container Requirements
- Docker requires root access (or rootless Docker)
- HPC systems typically use Singularity instead of Docker
- Singularity can convert Docker images. The toolkit images are built locally from each skill's `scripts/Dockerfile`, so convert from the local daemon and pass the result with `--container`:
  `singularity build pipeline-chipseq.sif docker-daemon://encode-toolkit/pipeline-chipseq:1.0.0`

### Genome Index Files
- STAR genome index requires ~32 GB RAM to generate and ~30 GB disk
- BWA index is smaller (~8 GB for human genome)
- Pre-built indices are available from ENCODE or iGenomes
- Log the exact index version and source in provenance

### Cloud Costs
- Forgot to stop instances = runaway costs
- Use preemptible/spot instances for 60–80% cost savings (with retry logic)
- Set billing alerts before starting cloud runs

### Resume and Checkpointing
- Always use `-resume` flag with Nextflow to avoid re-running completed steps
- Cromwell provides similar call caching
- This is critical for long-running pipelines (WGBS, Hi-C)

## Child Pipeline Skills

For detailed, executable pipeline implementations, use these assay-specific child skills:

| Pipeline Skill | Assay | Aligner | Caller |
|---------------|-------|---------|--------|
| `pipeline-chipseq` | ChIP-seq | BWA-MEM | MACS2 + IDR |
| `pipeline-atacseq` | ATAC-seq | Bowtie2 | MACS2 (Tn5-adjusted) |
| `pipeline-rnaseq` | RNA-seq | STAR | RSEM + Kallisto |
| `pipeline-wgbs` | WGBS | Bismark | MethylDackel |
| `pipeline-hic` | Hi-C | BWA | Juicer + HiCCUPS |
| `pipeline-dnaseseq` | DNase-seq | BWA | Hotspot2 |
| `pipeline-cutandrun` | CUT&RUN | Bowtie2 | SEACR |

Each child includes: SKILL.md overview, 5 stage reference files, Nextflow DSL2 pipeline, Dockerfile, and cloud deployment configs (local/SLURM/GCP/AWS).

## Walkthrough: Selecting and Configuring the Right Pipeline for Your ENCODE Data

**Goal**: Guide a researcher from raw ENCODE FASTQ files through pipeline selection, configuration, and execution using the appropriate ENCODE uniform processing pipeline.
**Context**: ENCODE provides standardized pipelines for each assay type. This skill helps users select the right pipeline and configure it for their specific experiment.

### Step 1: Identify the experiment and assay type

```
encode_get_experiment(accession="ENCSR000AKA")
```

Expected output:
```json
{
  "accession": "ENCSR000AKA",
  "assay_title": "Histone ChIP-seq",
  "target": "H3K27ac",
  "biosample_summary": "GM12878",
  "bio_replicate_count": 2,
  "assembly": ["GRCh38"],
  "status": "released"
}
```

**Interpretation**: This is a Histone ChIP-seq experiment targeting H3K27ac. Use the **pipeline-chipseq** skill for processing.

### Step 2: Download raw FASTQ files

```
encode_list_files(experiment_accession="ENCSR000AKA", file_format="fastq")
```

Expected output (a JSON array of files; fields abridged):
```json
[
  {"accession": "ENCFF001FQ1", "output_type": "reads", "file_format": "fastq", "biological_replicates": [1], "file_size_human": "2.3 GB", "status": "released"},
  {"accession": "ENCFF002FQ2", "output_type": "reads", "file_format": "fastq", "biological_replicates": [1], "file_size_human": "2.4 GB", "status": "released"}
]
```

### Step 3: Select pipeline based on assay type

| ENCODE Assay | Pipeline Skill | Key Tool |
|---|---|---|
| Histone ChIP-seq | **pipeline-chipseq** | BWA-MEM + MACS2 + IDR |
| TF ChIP-seq | **pipeline-chipseq** | BWA-MEM + MACS2 + IDR |
| ATAC-seq | **pipeline-atacseq** | Bowtie2 + Tn5 shift + MACS2 |
| RNA-seq | **pipeline-rnaseq** | STAR 2-pass + RSEM |
| WGBS | **pipeline-wgbs** | Bismark + MethylDackel |
| Hi-C | **pipeline-hic** | BWA + pairtools + Juicer |
| DNase-seq | **pipeline-dnaseseq** | BWA + Hotspot2 |
| CUT&RUN/CUT&Tag | **pipeline-cutandrun** | Bowtie2 + SEACR |

### Step 4: Configure and run

For Histone ChIP-seq. `--reads` is a glob that Nextflow's `fromFilePairs` must resolve, so give
the downloaded ENCFF files `_R1`/`_R2` names first. Which mate each accession is comes from its
page on encodeproject.org (`paired_end` 1 or 2, and `paired_with` naming the other accession):

```bash
mkdir -p fastq
ln -s "$PWD/ENCFF001FQ1.fastq.gz" fastq/rep1_R1.fastq.gz
ln -s "$PWD/ENCFF002FQ2.fastq.gz" fastq/rep1_R2.fastq.gz

nextflow run pipeline-chipseq/scripts/main.nf \
  -profile local \
  --reads 'fastq/*_R{1,2}.fastq.gz' \
  --genome GRCh38 \
  --peak_type narrow \
  --bwa_index ./GRCh38_index \
  --chrom_sizes /ref/hg38.chrom.sizes \
  --outdir results/ \
  -resume
```

`--reads` and `--chrom_sizes` are required. There is no `--target` parameter: the ChIP target
does not change the workflow, only `--peak_type` (`narrow` or `broad`, one value per run) does.
Add `--control '<glob>'` for input samples, making sure the two globs do not match the same
files.

### Step 5: Quality check the output

Use → **quality-assessment** skill to evaluate pipeline output against ENCODE standards:
- FRiP >= 1%
- NSC > 1.05
- RSC > 0.8

`pipeline-chipseq` computes FRiP for every treatment sample and publishes it as
`qc/<sample>.frip_mqc.tsv` (also a MultiQC table). NSC and RSC are not computed:
phantompeakqualtools is not in the pipeline image, so those stay manual steps on the filtered
BAM. The MultiQC report covers FastQC, trimming, flagstat, duplication metrics and FRiP.

### Integration with downstream skills
- Raw data from → **download-encode** provides FASTQ input for all pipelines
- Pipeline output feeds into → **quality-assessment** for ENCODE-standard QC
- Processed peaks feed into → **peak-annotation**, **regulatory-elements**, **histone-aggregation**
- Each assay has a dedicated pipeline skill: pipeline-chipseq through pipeline-cutandrun

## Code Examples

### 1. Determine which pipeline to use
```
encode_get_experiment(accession="ENCSR000AKA")
```

Expected output:
```json
{
  "accession": "ENCSR000AKA",
  "assay_title": "Histone ChIP-seq",
  "assembly": ["GRCh38"],
  "target": "H3K27ac"
}
```

### 2. Find FASTQ files for pipeline input
```
encode_list_files(experiment_accession="ENCSR000AKA", file_format="fastq")
```

Expected output (a JSON array of files; fields abridged):
```json
[
  {"accession": "ENCFF001FQ1", "output_type": "reads", "file_format": "fastq", "file_size_human": "2.3 GB", "status": "released"},
  {"accession": "ENCFF002FQ2", "output_type": "reads", "file_format": "fastq", "file_size_human": "2.4 GB", "status": "released"}
]
```

### 3. Survey available data by assay type for pipeline selection
```
encode_get_facets(organism="Homo sapiens")
```

Expected output (top-level keys are ENCODE facet field names):
```json
{
  "assay_title": [
    {"term": "Histone ChIP-seq", "count": 2500},
    {"term": "TF ChIP-seq", "count": 1800},
    {"term": "total RNA-seq", "count": 1200},
    {"term": "ATAC-seq", "count": 450},
    {"term": "WGBS", "count": 147}
  ]
}
```

## Integration

| This skill produces... | Feed into... | Purpose |
|---|---|---|
| Pipeline selection recommendation | **pipeline-chipseq** through **pipeline-cutandrun** | Route to correct assay-specific pipeline |
| FASTQ download commands | **download-encode** | Obtain raw data for pipeline input |
| Pipeline configuration | **bioinformatics-installer** | Install required pipeline dependencies |
| Pipeline output files | **quality-assessment** | Validate output against ENCODE QC standards |
| Processed peaks/signals | **peak-annotation** | Annotate pipeline output with gene assignments |
| Processed peaks | **regulatory-elements** | Classify pipeline output as enhancers/promoters/insulators |
| Pipeline run metadata | **data-provenance** | Log pipeline parameters and versions |
| Processed data | **visualization-workflow** | Generate QC and analysis visualizations |

## Related Skills

- `data-provenance` — Exact provenance logging for every operation
- `quality-assessment` — Evaluating pipeline output quality
- `download-encode` — Downloading ENCODE files for pipeline input
- `single-cell-encode` — Single-cell pipeline specifics
- `publication-trust` — Verify literature claims backing analytical decisions

## Presenting Results

When reporting pipeline recommendations:

- **Selected pipeline**: State the recommended pipeline (e.g., pipeline-chipseq, pipeline-atacseq) with a brief rationale based on the assay type and user's data
- **Resource estimates**: Present CPU, RAM, disk, and estimated runtime requirements in a table, compared against the user's available resources from system checks
- **Container availability**: Confirm whether Docker or Singularity is available and report the recommended container image with version tag
- **Execution profile**: Recommend the appropriate profile (local, slurm, gcp, aws) based on the user's compute environment
- **Cost estimate**: For cloud execution, provide per-sample cost estimates and recommend preemptible/spot instances where applicable
- **Genome index status**: Note whether pre-built genome indices are available or need to be generated, and estimate the index build time
- **Configuration summary**: Provide the recommended nextflow.config parameters tailored to the user's system
- **Next steps**: Direct the user to the specific child pipeline skill (e.g., "Use `pipeline-chipseq` to execute the pipeline with the parameters above")

## For the request: "$ARGUMENTS"
