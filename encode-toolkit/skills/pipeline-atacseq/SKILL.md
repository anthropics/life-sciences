---
name: pipeline-atacseq
description: "Execute ENCODE ATAC-seq processing pipeline from FASTQ to peaks and signal tracks. Child of pipeline-guide. Provides stage-by-stage Nextflow execution with Docker containers and cloud deployment. Handles Tn5 transposase offset correction, mitochondrial read removal, and nucleosome-free fragment selection. Use when users need to process ATAC-seq data following ENCODE standards. Trigger on: ATAC-seq pipeline, run ATAC-seq, process ATAC-seq, chromatin accessibility, open chromatin, Tn5 shift, TSS enrichment."
---

# ENCODE ATAC-seq Pipeline

## When to Use

- User wants to run an ATAC-seq processing pipeline from FASTQ to peaks and signal tracks
- User asks about "ATAC-seq pipeline", "Tn5 shift", "chromatin accessibility pipeline", or "Bowtie2 for ATAC"
- User needs to process ATAC-seq data with proper Tn5 insertion site correction
- Example queries: "process my ATAC-seq FASTQs", "run ENCODE ATAC-seq pipeline", "call accessibility peaks from ATAC-seq"

Execute the ENCODE ATAC-seq processing pipeline from raw FASTQ files through Tn5 offset
correction, peak calling, IDR analysis, and signal track generation. This skill provides
a Nextflow DSL2 implementation following ENCODE uniform analysis standards.

TSS enrichment scoring is a **manual post-processing step**; the workflow does not compute
it (see "Manual QC steps" below and `references/05-qc-metrics.md`).

## Overview

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) uses the Tn5
transposase to probe open chromatin regions. This pipeline processes ATAC-seq data
through quality control, alignment with Bowtie2, mitochondrial read removal, duplicate
removal, Tn5 insertion site correction (+4/-5 bp offset), blacklist filtering,
nucleosome-free fragment selection, MACS2 peak calling, FRiP calculation, and an IDR
comparison for every pair of replicates.

Key differences from ChIP-seq: Bowtie2 aligner (optimized for short fragments), Tn5
transposase shift correction, mitochondrial read filtering (chrM can be 30-80% of reads),
and no input control.

The workflow is **paired-end only**. Passing `--single_end` stops the run with an error,
because Tn5 shifting, nucleosome-free selection and BAMPE peak calling all depend on
fragment length.

## Key Literature

| Reference | Journal | Year | DOI | Relevance |
|-----------|---------|------|-----|-----------|
| Buenrostro et al. "Transposition of native chromatin (ATAC-seq)" | Nature Methods | 2013 | 10.1038/nmeth.2688 | Original ATAC-seq method (~5,000 citations) |
| Corces et al. "An improved ATAC-seq protocol" | Nature Methods | 2017 | 10.1038/nmeth.4396 | Omni-ATAC improvements (~2,500 citations) |
| ENCODE Project Consortium "Expanded encyclopaedias" | Nature | 2020 | 10.1038/s41586-020-2493-4 | ENCODE Phase 3 standards |
| Amemiya et al. "ENCODE Blacklist" | Scientific Reports | 2019 | 10.1038/s41598-019-45839-z | Artifact regions (~1,372 citations) |
| Langmead & Salzberg "Fast gapped-read alignment with Bowtie 2" | Nature Methods | 2012 | 10.1038/nmeth.1923 | Aligner (~30,000 citations) |
| Yan et al. "From reads to insight: ATAC-seq analysis" | Genome Biology | 2020 | 10.1186/s13059-020-1929-3 | Analysis best practices |

## Pipeline Stages

```
FASTQ ──> FastQC / Trim Galore ──> Bowtie2 ──> Mito Removal ──> Picard MarkDuplicates
  │                                            (chrM dropped)   (duplicates REMOVED)
  │                                                                       │
  │           ┌───────────────────────────────────────────────────────────┘
  │           v
  │     Tn5 Shift (alignmentSieve --ATACshift) ──> Blacklist Filter ──> Size Selection
  │                                                       │                    │
  │                                                       │        ┌───────────┴────────┐
  │                                                       v        v                    v
  │                                               Signal Track   NFR (<150 bp)   Mono-nucleosome
  │                                                (all frags)      │              (150-300 bp)
  │                                                                 v
  │                                                   MACS2 Peak Calling ──> IDR (every pair)
  │                                                                 │
  │                                                                 v
  │                                                   FRiP (NFR peaks vs final BAM)
  v
 QC reports ────────────────────────────────────────────────────────────────> MultiQC
```

The Tn5 shift runs **after** duplicate removal, and peaks are called on the
nucleosome-free BAM only. The signal track is built from all fragments in the
blacklist-filtered BAM, not from the NFR BAM.

### Stage Summary

| Stage | Tool | Input | Output | Reference |
|-------|------|-------|--------|-----------|
| 1. QC & Trimming | FastQC, Trim Galore | Raw FASTQ | Trimmed FASTQ, FastQC reports | references/01-qc-trimming.md |
| 2. Alignment | Bowtie2, samtools | Trimmed FASTQ | Sorted BAM, flagstat, bowtie2 log | references/02-alignment.md |
| 3. Filtering & Tn5 shift | samtools, Picard, deeptools `alignmentSieve`, bedtools | Sorted BAM | Shifted, filtered, size-selected BAMs | references/03-tn5-filtering.md |
| 4. Peak Calling & IDR | MACS2, IDR | NFR BAM | narrowPeak, one `<sampleA>_vs_<sampleB>.idr_peaks.txt` per replicate pair | references/04-peak-calling.md |
| 5. Signal, FRiP & QC report | deeptools `bamCoverage`, bedtools, samtools, MultiQC | Filtered BAM, NFR peaks, QC logs | bigWig, `<sample>.frip_mqc.tsv`, multiqc_report.html | references/05-qc-metrics.md |

## Input Requirements

### Required
- **ATAC-seq FASTQ** (`--reads`): paired-end reads, gzipped. A Nextflow file-pair glob,
  e.g. `'fastq/*_R{1,2}.fq.gz'`.
- **Bowtie2 index directory** (`--bowtie2_index`): Bowtie2 is invoked as
  `bowtie2 ... -x <dir>/<genome>`, so the directory must hold index files named after the
  genome:

```
GRCh38_bowtie2_index/
  GRCh38.1.bt2  GRCh38.2.bt2  GRCh38.3.bt2  GRCh38.4.bt2
  GRCh38.rev.1.bt2  GRCh38.rev.2.bt2
```

  Build it once with `bowtie2-build GRCh38.fa GRCh38_bowtie2_index/GRCh38`. If the flag is
  omitted, the workflow looks for `./<genome>_bowtie2_index` in the launch directory. The
  workflow does not build or download the index.

### Optional
- **Blacklist** (`--blacklist`): defaults to the ENCODE Blacklist v2 URL for `--genome`.

There is no sample sheet and no input control. Inputs are globs, and every sample in a run
shares one `--genome`. Unlike ChIP-seq, ATAC-seq does not need a separate input or IgG
control; MACS2 calls peaks against a local background model.

## Tn5 Transposase Offset Correction

The Tn5 transposase inserts sequencing adapters with a 9-bp duplication. To center
reads on the actual cut site:
- **Forward strand (+)**: shift +4 bp
- **Reverse strand (-)**: shift -5 bp

The workflow applies this with `alignmentSieve --ATACshift` (deeptools) after duplicate
removal and before blacklist filtering. The correction is essential for footprinting and
motif analysis.

## Fragment Size Distribution

ATAC-seq produces a characteristic nucleosomal ladder pattern:

| Fragment Class | Size Range | Biological Meaning |
|---------------|------------|-------------------|
| Nucleosome-free (NFR) | <150 bp | Open chromatin / TF binding |
| Mono-nucleosome | 150-300 bp | Single nucleosome wrapping |
| Di-nucleosome | 300-500 bp | Two nucleosomes |
| Tri-nucleosome | 500-700 bp | Three nucleosomes |

The workflow calls peaks on the nucleosome-free BAM. The NFR/mono-nucleosome boundary is
`--nfr_max` (default 150); the mono-nucleosome selection is `--nfr_max` to 300 bp. The
workflow does not plot the fragment size distribution.

## Parameters

### Pipeline parameters (`main.nf`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | none (required) | Glob for the paired-end FASTQ file pairs |
| `--bowtie2_index` | `./<genome>_bowtie2_index` | Directory holding the Bowtie2 index files named `<genome>.*.bt2` |
| `--genome` | `GRCh38` | `GRCh38` or `mm10`; sets the MACS2 genome size, the default index directory and the default blacklist |
| `--blacklist` | ENCODE Blacklist v2 URL for `--genome` | BED (or `.bed.gz`) of artifact regions removed from the BAM |
| `--mito_name` | `chrM` | Name of the mitochondrial contig to drop |
| `--nfr_max` | `150` | Maximum nucleosome-free fragment length, and the lower bound of the mono-nucleosome selection |
| `--skip_idr` | `false` | Skip the IDR step |
| `--single_end` | `false` | Accepted but always rejected: the workflow stops with an error because it is paired-end only |
| `--outdir` | `results` | Where results are published |

### Infrastructure parameters (`nextflow.config`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--container` | `encode-toolkit/pipeline-atacseq:1.0.0` | Image built from `scripts/Dockerfile`. Pass a registry image for `gcp`/`aws`, or a `.sif` file for `slurm` |
| `--max_cpus`, `--max_memory`, `--max_time` | `16`, `64.GB`, `24.h` | Upper bounds applied to every process |
| `--slurm_queue`, `--slurm_account` | `normal`, none | SLURM partition and account |
| `--gcp_project`, `--gcp_workdir` | none (both required for `-profile gcp`) | Google Cloud project and `gs://` work directory |
| `--gcp_location`, `--gcp_disk` | `us-central1`, `200.GB` | Google Batch region and per-task disk |
| `--aws_queue`, `--aws_workdir` | none (both required for `-profile aws`) | AWS Batch job queue and `s3://` work directory |
| `--aws_region`, `--aws_cli_path` | `us-east-1`, `/home/ec2-user/miniconda/bin/aws` | AWS region, and the AWS CLI path inside the Batch AMI |

Profiles are `local`, `slurm`, `gcp` and `aws`.

## QC Thresholds

**The workflow computes only the metrics marked "workflow" below.** TSS enrichment,
NRF/PBC, fragment-size plots and ataqv are manual post-processing steps documented in
`references/05-qc-metrics.md`.

| Metric | Threshold | Computed by | Source |
|--------|-----------|-------------|--------|
| Total sequenced reads | >=50M (recommended) | workflow (FastQC, flagstat) | ENCODE |
| Mapping rate | >80% | workflow (bowtie2 log, `samtools flagstat`) | ENCODE |
| Mitochondrial fraction | <20% (ideal <5%) | workflow (`qc/<sample>.idxstats.txt`) | ENCODE |
| Duplication rate | <30% | workflow (Picard `dup_metrics.txt`) | ENCODE |
| IDR peaks at 0.05 | >50,000 | workflow (`peaks/idr/<sampleA>_vs_<sampleB>.idr_peaks.txt`) | ENCODE |
| NRF (non-redundant fraction) | >=0.8 | manual | ENCODE |
| PBC1 | >=0.8 | manual | ENCODE |
| TSS enrichment score | >=5 (GRCh38), >=6 (hg19), >=10 (mm10) | manual (deeptools + a TSS BED) | ENCODE standard |
| FRiP | >=0.3 | workflow (`qc/<sample>.frip_mqc.tsv`) | ENCODE |
| NFR fraction | >0.4 of fragments <150bp | manual | Buenrostro 2013 |

`qc/<sample>.idxstats.txt` is `samtools idxstats` of the BAM before mitochondrial reads are
removed (contig, length, mapped, unmapped): the mitochondrial fraction is the mapped count
on the `--mito_name` row divided by the sum of the mapped column. MultiQC's samtools module
reads the same file and reports that fraction.

### TSS Enrichment Score (manual)

The TSS enrichment score measures the fold enrichment of ATAC-seq signal at
transcription start sites compared to flanking regions. It is the single most
informative QC metric for ATAC-seq, but **this workflow does not compute it**: there is no
TSS BED input and no `computeMatrix`/`plotProfile` step. Run it manually against
`signal/<sample>.signal.bw` with a TSS BED for your assembly; the commands are in
`references/05-qc-metrics.md`.

| Score | Quality | Interpretation |
|-------|---------|---------------|
| >=7 | Excellent | High signal-to-noise |
| 5-7 | Good | Acceptable for most analyses |
| 3-5 | Marginal | Review other metrics carefully |
| <3 | Poor | Likely failed; consider re-doing |

## Execution

### Quick Start (Local Docker)
```bash
nextflow run scripts/main.nf \
  -profile local \
  --reads 'fastq/*_R{1,2}.fq.gz' \
  --genome GRCh38 \
  --bowtie2_index GRCh38_bowtie2_index \
  --blacklist hg38-blacklist.v2.bed.gz \
  --outdir results/
```

`--blacklist` is optional; without it the workflow downloads the ENCODE Blacklist v2 for
`--genome`. Give the glob at least two replicates if you want the IDR step to run.

### SLURM HPC

The `slurm` profile runs through Singularity, so pass a local image file rather than the
default Docker image name:

```bash
singularity build pipeline-atacseq.sif docker-daemon://encode-toolkit/pipeline-atacseq:1.0.0

nextflow run scripts/main.nf \
  -profile slurm \
  --container /path/to/pipeline-atacseq.sif \
  --slurm_queue normal \
  --reads 'fastq/*_R{1,2}.fq.gz' \
  --genome GRCh38 \
  --bowtie2_index GRCh38_bowtie2_index \
  --outdir results/
```

### Cloud

```bash
# Google Cloud Batch
nextflow run scripts/main.nf -profile gcp \
    --container us-docker.pkg.dev/<project>/<repo>/pipeline-atacseq:1.0.0 \
    --gcp_project <project> \
    --gcp_workdir gs://<bucket>/work \
    --reads 'gs://<bucket>/fastq/*_R{1,2}.fq.gz' \
    --genome GRCh38 \
    --bowtie2_index gs://<bucket>/reference/GRCh38_bowtie2_index \
    --outdir gs://<bucket>/results

# AWS Batch
nextflow run scripts/main.nf -profile aws \
    --container <account>.dkr.ecr.<region>.amazonaws.com/pipeline-atacseq:1.0.0 \
    --aws_queue <job-queue> \
    --aws_workdir s3://<bucket>/work \
    --reads 's3://<bucket>/fastq/*_R{1,2}.fq.gz' \
    --genome GRCh38 \
    --bowtie2_index s3://<bucket>/reference/GRCh38_bowtie2_index \
    --outdir s3://<bucket>/results
```

`--outdir` only sets where results are published; Google Batch and AWS Batch stage every
task through the work directory, and the workflow stops with an error if it or the
project/queue is missing.

## Cloud Cost Estimates

| Platform | Instance | Cost/Sample | Time/Sample | Notes |
|----------|----------|-------------|-------------|-------|
| GCP | n1-standard-8 | ~$2-4 | 2-3 hours | Spot VMs enabled in the `gcp` profile |
| AWS | m5.2xlarge | ~$2-4 | 2-3 hours | Spot instances recommended |
| Local | 8 cores, 32GB | $0 | 3-5 hours | Docker required |
| SLURM | 8 cores, 32GB | Varies | 2-3 hours | Singularity image required |

## Output Directory Structure

```
results/
  fastqc/                   # FastQC reports for raw and trimmed reads (.html, .zip)
  trimmed/                  # Trimmed FASTQ (*_val_1.fq.gz / *_val_2.fq.gz) + trimming reports
  aligned/                  # <sample>.bam, .bam.bai, <sample>.flagstat.txt,
                            #   <sample>.bowtie2.log
  filtered/                 # <sample>.dup_metrics.txt, <sample>.final.bam(.bai),
                            #   <sample>.final.flagstat.txt
    shifted/                # <sample>.shifted.bam(.bai) -- Tn5-corrected, pre-blacklist
    nfr/                    # <sample>.nfr.bam(.bai) and <sample>.mononuc.bam
  peaks/
    narrow/                 # <sample>_peaks.narrowPeak, _summits.bed, _peaks.xls,
                            #   _treat_pileup.bdg, _control_lambda.bdg
    idr/                    # <sampleA>_vs_<sampleB>.idr_peaks.txt (+ .png), one file per
                            #   replicate pair
  signal/                   # <sample>.signal.bw (all fragments, RPKM)
  qc/
    <sample>.idxstats.txt   # samtools idxstats before chrM removal:
                            #   contig / length / mapped / unmapped
    <sample>.frip_mqc.tsv   # Peak set / FRiP / reads_in_peaks / total_reads
    multiqc/                # multiqc_report.html, multiqc_data/
  pipeline_info/            # timeline.html, report.html, trace.txt
```

The nucleosome-free and mono-nucleosome BAMs are both written to `filtered/nfr/`; there is
no `filtered/mononuc/` directory, and `mononuc.bam` has no index.

## Common Pitfalls

### 1. High Mitochondrial Read Fraction
Mitochondrial DNA lacks chromatin and is highly accessible, often capturing 30-80%
of reads. This is the most common ATAC-seq quality issue. The workflow removes `--mito_name`
reads and publishes the per-contig counts it used, `qc/<sample>.idxstats.txt`: divide the
mapped count on the `--mito_name` row by the sum of the mapped column to get the fraction,
or read it from the samtools section of the MultiQC report. If >50% mito, consider
optimizing the cell lysis step.

### 2. Wrong mitochondrial contig name
`--mito_name` defaults to `chrM`. Assemblies that call the contig `MT` need
`--mito_name MT`. Confirm the name with `samtools idxstats` on a BAM from your index
before running.

### 3. Using BWA Instead of Bowtie2
Bowtie2 handles the short fragments from ATAC-seq (especially NFR <150bp) better
than BWA-MEM. The workflow uses Bowtie2 with `--very-sensitive`.

### 4. Only one replicate in the `--reads` glob
IDR needs a pair of samples. With one sample there is no pair, so IDR is skipped silently
and `peaks/idr/` is never created. Every sample matched by `--reads` is treated as a
replicate of the same experiment, so unrelated samples in one glob produce meaningless
pairwise comparisons.

### 5. TSS enrichment is not in the output
TSS enrichment is the most informative single metric for ATAC-seq, but the workflow does
not compute it. Run the manual `computeMatrix`/`plotProfile` step in
`references/05-qc-metrics.md` before judging a library.

## Pipeline Scripts

| File | Description |
|------|-------------|
| `scripts/main.nf` | Nextflow DSL2 pipeline |
| `scripts/nextflow.config` | Execution profiles (local/slurm/gcp/aws) |
| `scripts/Dockerfile` | Docker image with all pipeline tools |

The image is pinned to `linux/amd64`; on an arm64 host it runs under emulation.

Tool versions in the image: Bowtie2 2.5.4, samtools 1.19, bedtools 2.31.0, Picard 3.1.1
(Java 17), Trim Galore 0.6.10 with cutadapt 4.6, FastQC 0.12.1, MACS2 2.2.9.1,
IDR 2.0.4.2, deepTools 3.5.5, MultiQC 1.21. The conda environment
`bioinformatics-installer/environments/atacseq-env.yml` pins the same versions of the
tools it lists.

## ENCODE Data Integration

After running on your own data, compare with ENCODE reference:

```python
# Find matching ENCODE ATAC-seq experiments
encode_search_experiments(
    assay_title="ATAC-seq",
    organ="pancreas",
    biosample_type="tissue"
)

# Download ENCODE peaks for comparison
encode_batch_download(
    download_dir="/data/encode_reference/",
    output_type="IDR thresholded peaks",
    assay_title="ATAC-seq",
    organ="pancreas",
    assembly="GRCh38"
)
```

## Pitfalls & Edge Cases

- **Tn5 shift is critical**: ATAC-seq reads must be shifted +4/-5 bp to center on the Tn5
  insertion site. The workflow does this with `alignmentSieve --ATACshift` after duplicate
  removal. Without the correction, footprinting is offset by ~5 bp.
- **Mitochondrial reads dominate**: Expect 30-80% mitochondrial reads. The workflow drops
  them right after alignment, before duplicate removal and peak calling. >80% chrM
  indicates dead/dying cells or poor nuclei isolation.
- **Fragment size distribution is diagnostic**: a nucleosomal ladder (sub-nucleosomal
  <150bp, mono-nucleosomal ~200bp, di-nucleosomal ~400bp) confirms successful
  transposition. The workflow does not plot it; use `bamPEFragmentSize` manually.
- **TSS enrichment threshold**: ENCODE requires TSS enrichment >=5 (GRCh38), >=6 (hg19), or
  >=10 (mm10) for ATAC-seq (ENCODE data standards). Values below 4 indicate poor
  signal-to-noise. Computed manually, not by this workflow.
- **MACS2 `--shift`/`--extsize` do not apply here**: the workflow calls peaks in `-f BAMPE`
  mode, where MACS2 takes fragment coordinates from read pairs, forces `--nomodel` and
  neutralises `--shift` internally. Shift/extension values only matter when calling peaks
  on BED or single-end input. See `references/04-peak-calling.md`.
- **Paired-end only**: `--single_end` is rejected with an error. Single-end ATAC-seq cannot
  distinguish nucleosome-free from nucleosomal fragments.

## Walkthrough: Processing ENCODE ATAC-seq from FASTQ to Accessible Chromatin Peaks

**Goal**: Process raw ATAC-seq FASTQ files through this pipeline to generate
nucleosome-free region peaks and a signal track.
**Context**: Bowtie2 alignment, chrM removal, duplicate removal, Tn5 shift (+4/-5),
blacklist filtering, NFR selection and MACS2 peak calling.

### Step 1: Find ATAC-seq experiment

```
encode_get_experiment(accession="ENCSR637ENO")
```

Expected output (fields abridged):
```json
{
  "accession": "ENCSR637ENO",
  "assay_title": "ATAC-seq",
  "biosample_summary": "GM12878",
  "assembly": ["GRCh38"],
  "bio_replicate_count": 2,
  "tech_replicate_count": 2,
  "status": "released"
}
```

### Step 2: List FASTQ files

```
encode_list_files(experiment_accession="ENCSR637ENO", file_format="fastq")
```

Expected output (a JSON array of files; fields abridged):
```json
[
  {"accession": "ENCFF100ATQ", "file_format": "fastq", "output_type": "reads", "biological_replicates": [1], "file_size_human": "1.8 GB"},
  {"accession": "ENCFF101ATQ", "file_format": "fastq", "output_type": "reads", "biological_replicates": [1], "file_size_human": "1.9 GB"},
  {"accession": "ENCFF102ATQ", "file_format": "fastq", "output_type": "reads", "biological_replicates": [2], "file_size_human": "1.7 GB"},
  {"accession": "ENCFF103ATQ", "file_format": "fastq", "output_type": "reads", "biological_replicates": [2], "file_size_human": "1.8 GB"}
]
```

The listing does not say which file of a pair is read 1 and which is read 2 -- no
`encode_*` tool reports that. Open each file's page on encodeproject.org, where
`paired_end` is 1 or 2 and `paired_with` names the other accession. Both replicates are
needed for the IDR step.

### Step 3: Download and name the FASTQs so a read-pair glob can find them

```
encode_download_files(file_accessions=["ENCFF100ATQ", "ENCFF101ATQ", "ENCFF102ATQ", "ENCFF103ATQ"], download_dir="/data/atacseq/fastq")
```

ENCODE FASTQs are named by accession (`ENCFF123ABC.fastq.gz`) with no `_R1`/`_R2` in the
name, so the `--reads` glob (`*_R{1,2}.fq.gz`) cannot pair them. Take the mate assignment
from each file's page on encodeproject.org (`paired_end` is 1 or 2, `paired_with` names the
other accession), then link them into the shape the glob expects:

```bash
cd /data/atacseq/fastq
ln -s ENCFF100ATQ.fastq.gz gm12878_rep1_R1.fq.gz
ln -s ENCFF101ATQ.fastq.gz gm12878_rep1_R2.fq.gz
ln -s ENCFF102ATQ.fastq.gz gm12878_rep2_R1.fq.gz
ln -s ENCFF103ATQ.fastq.gz gm12878_rep2_R2.fq.gz
```

### Step 4: Run the ATAC-seq pipeline

```bash
nextflow run scripts/main.nf \
  -profile local \
  --reads '/data/atacseq/fastq/gm12878_*_R{1,2}.fq.gz' \
  --genome GRCh38 \
  --bowtie2_index /data/reference/GRCh38_bowtie2_index \
  --blacklist /data/reference/hg38-blacklist.v2.bed.gz \
  --mito_name chrM \
  --outdir /data/atacseq/results
```

Pipeline steps, in the order the workflow runs them:
1. FastQC on raw reads
2. Adapter trimming with Trim Galore (`--nextera`), plus FastQC on the trimmed reads
3. Alignment (Bowtie2 `--very-sensitive`, MAPQ 30, properly paired only)
4. Mitochondrial read removal, with `samtools idxstats` of the pre-removal BAM published as `qc/<sample>.idxstats.txt`
5. Duplicate removal (Picard `REMOVE_DUPLICATES=true`)
6. Tn5 shift correction (+4/-5, `alignmentSieve --ATACshift`)
7. Blacklist filtering of the BAM
8. Nucleosome-free (<150 bp) and mono-nucleosome (150-300 bp) selection
9. Peak calling on the NFR BAM (MACS2 `-f BAMPE --nomodel --keep-dup all --call-summits --qvalue 0.05 -B`)
10. IDR on every pair of replicates, signal track from all fragments, FRiP, MultiQC

### Step 5: Validate output quality

From the workflow:
| Output | What to check |
|---|---|
| `qc/multiqc/multiqc_report.html` | Mapping rate (>80%), adapter content, duplication rate, mitochondrial fraction, FRiP |
| `qc/<sample>.idxstats.txt` | Mitochondrial fraction (<20%, ideal <5%): mapped reads on the `chrM` row over the sum of the mapped column |
| `qc/<sample>.frip_mqc.tsv` | FRiP (>=0.3 for ATAC-seq) |
| `peaks/narrow/<sample>_peaks.narrowPeak` | Peak count per replicate |
| `peaks/idr/gm12878_rep1_vs_gm12878_rep2.idr_peaks.txt` | IDR peaks at 0.05 (>50,000) |

With the two replicates above there is one IDR file; a third replicate would add
`gm12878_rep1_vs_gm12878_rep3` and `gm12878_rep2_vs_gm12878_rep3`.

Manual follow-ups (not run by this workflow): TSS enrichment, fragment-size distribution
plots, NRF/PBC and ataqv. Commands are in `references/05-qc-metrics.md`.

### Step 6: Track and log provenance

```
encode_track_experiment(accession="ENCSR637ENO", notes="GM12878 ATAC-seq processed through the pipeline-atacseq skill")
```

### Integration with downstream skills
- Accessible chromatin peaks feed into -> **accessibility-aggregation** for cross-experiment union merge
- Peak regions feed into -> **motif-analysis** for TF motif enrichment
- Signal tracks feed into -> **visualization-workflow** for browser display
- Peaks feed into -> **regulatory-elements** for cCRE classification
- QC metrics validated by -> **quality-assessment**

## Code Examples

### 1. Find ATAC-seq data for processing

```
encode_search_experiments(
  assay_title="ATAC-seq",
  organ="pancreas"
)
```

Expected output:
```json
{
  "results": [
    {
      "accession": "ENCSR789PAN",
      "assay_title": "ATAC-seq",
      "biosample_summary": "pancreas tissue male adult (44 years)",
      "status": "released"
    }
  ],
  "total": 8,
  "limit": 25,
  "offset": 0,
  "has_more": false,
  "next_offset": null
}
```

### 2. Check file details before download

```
encode_list_files(
  experiment_accession="ENCSR789PAN",
  file_format="fastq"
)
```

Expected output:
```json
[
  {
    "accession": "ENCFF100ATQ",
    "file_format": "fastq",
    "output_type": "reads",
    "biological_replicates": [1],
    "file_size_human": "3.1 GB",
    "status": "released"
  }
]
```

## Integration

| This skill produces... | Feed into... | Purpose |
|---|---|---|
| Accessible chromatin peaks | **accessibility-aggregation** | Cross-experiment union merge |
| Peak regions (BED) | **motif-analysis** | TF motif enrichment in open chromatin |
| Signal tracks (bigWig) | **visualization-workflow** | Genome browser accessibility display |
| Nucleosome-free peaks | **regulatory-elements** | Classify accessible regions as enhancers/promoters |
| Peak coordinates | **variant-annotation** | Identify variants in accessible chromatin |
| QC outputs (`idxstats.txt`, `frip_mqc.tsv`, MultiQC) | **quality-assessment** | Validate against ENCODE ATAC-seq standards |
| `pipeline_info/` reports | **data-provenance** | Record Tn5 shift, fragment filters, tool versions |
| Peak files | **jaspar-motifs** | Scan accessible regions for known TF motifs |

## Related Skills

- **pipeline-guide** (parent): General pipeline selection and resource assessment
- **accessibility-aggregation**: Merge ATAC-seq peaks across samples
- **quality-assessment**: Deep-dive QC analysis beyond basic metrics
- **regulatory-elements**: Annotate peaks with regulatory element classifications
- **compare-biosamples**: Compare accessibility profiles across cell types
- **pipeline-chipseq**: Sibling pipeline for ChIP-seq data
- **publication-trust**: Verify literature claims backing analytical decisions

## Presenting Results

When reporting ATAC-seq pipeline results:

- **Mitochondrial fraction**: Report it from `qc/<sample>.idxstats.txt` -- mapped reads on
  the `--mito_name` row over the sum of the mapped column -- or from the samtools section
  of the MultiQC report (ideal <5%, acceptable <20%)
- **Key QC metrics from the run**: mapping rate (bowtie2 log, `samtools flagstat`),
  duplication rate (Picard `dup_metrics.txt`), and read counts, all aggregated in
  `qc/multiqc/multiqc_report.html`
- **Peak counts**: Report the per-replicate MACS2 peak count and, for every replicate pair,
  the IDR peak count at the 0.05 threshold. IDR is skipped when only one sample was
  processed
- **TSS enrichment**: State plainly that the workflow does not compute it. Report it only
  if the user ran the manual step, with the quality tier (Excellent >=7, Good 5-7,
  Marginal 3-5, Poor <3)
- **FRiP**: Report the value from `qc/<sample>.frip_mqc.tsv` (>=0.3 for ATAC-seq). It is
  computed from the blacklist-filtered BAM (all fragments) against the peaks called on the
  nucleosome-free fragments
- **Fragment size distribution / NFR fraction / ataqv**: also manual; do not report values
  the run did not produce
- **Output paths**: List key outputs (`peaks/narrow/`, `peaks/idr/`, `signal/`,
  `filtered/nfr/`, `qc/`, `pipeline_info/`)
- **Next steps**: Suggest `motif-analysis` for TF footprinting and de novo motif discovery,
  or `visualization-workflow` for genome browser session generation

## For the request: "$ARGUMENTS"
