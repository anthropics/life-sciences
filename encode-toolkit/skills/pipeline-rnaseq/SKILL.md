---
name: pipeline-rnaseq
description: "Execute ENCODE RNA-seq pipeline from FASTQ to gene quantification and signal tracks. Child of pipeline-guide. Provides Nextflow execution with Docker and cloud deployment. Use when processing RNA-seq data with STAR alignment, RSEM/Kallisto quantification, or generating expression matrices. Trigger on: RNA-seq pipeline, gene expression, STAR alignment, RSEM quantification, transcript quantification, TPM, FPKM, RNA processing, run RNA-seq."
---

# ENCODE RNA-seq Pipeline

## When to Use

- User wants to run an RNA-seq processing pipeline from FASTQ to gene quantification
- User asks about "RNA-seq pipeline", "STAR alignment", "RSEM", "gene expression quantification", or "Kallisto"
- User needs to process bulk RNA-seq data with ENCODE-standard 2-pass STAR alignment
- Example queries: "process my RNA-seq FASTQs", "quantify gene expression from RNA-seq", "run STAR and RSEM on my data"

Execute the ENCODE RNA-seq processing pipeline from raw FASTQ files through splice-aware
alignment, gene/transcript quantification, and strand-specific signal track generation.
This skill provides a complete Nextflow DSL2 implementation following ENCODE uniform
analysis standards.

## Overview

RNA-seq measures transcriptome-wide gene expression by sequencing cDNA derived from
cellular RNA. The ENCODE pipeline processes RNA-seq data through quality control,
splice-aware alignment with STAR (2-pass mode), gene and transcript quantification
with RSEM, optional fast pseudoalignment with Kallisto, and generation of strand-specific
signal tracks as bigWig files.

Key design decisions: STAR 2-pass mode for maximum splice junction sensitivity, RSEM
for accurate gene/transcript/isoform quantification including multi-mapped reads,
stranded library protocol (dUTP/rf-stranded) as the ENCODE standard, and paired-end
sequencing with a minimum of 30 million uniquely mapped reads per replicate.

## Key Literature

| Reference | Journal | Year | DOI | Relevance |
|-----------|---------|------|-----|-----------|
| Dobin et al. "STAR: ultrafast universal RNA-seq aligner" | Bioinformatics | 2013 | 10.1093/bioinformatics/bts635 | Splice-aware aligner (~12,000 citations) |
| Li & Dewey "RSEM: accurate transcript quantification from RNA-Seq data" | BMC Bioinformatics | 2011 | 10.1186/1471-2105-12-323 | Gene/transcript quantification (~6,000 citations) |
| Bray et al. "Near-optimal probabilistic RNA-seq quantification" | Nature Biotechnology | 2016 | 10.1038/nbt.3519 | Fast pseudoalignment (~4,000 citations) |
| Wang et al. "RSeQC: quality control of RNA-seq experiments" | Bioinformatics | 2012 | 10.1093/bioinformatics/bts356 | RNA-seq QC suite (~3,500 citations) |
| ENCODE Project Consortium "Expanded encyclopaedias" | Nature | 2020 | 10.1038/s41586-020-2493-4 | ENCODE Phase 3 standards |
| Frankish et al. "GENCODE 2021" | Nucleic Acids Research | 2021 | 10.1093/nar/gkaa1087 | Gene annotation reference |

## Pipeline Stages

```
FASTQ
  ├─> FastQC (raw reads)
  └─> Trim Galore (+ FastQC on the trimmed reads)
        ├─> Kallisto (optional) ──────────> kallisto/<sample>/abundance.tsv
        └─> STAR (2-pass)
              ├─> transcriptome BAM ─> RSEM ─> <sample>.genes.results / .isoforms.results
              ├─> bedGraph str1/str2 ─> bedGraphToBigWig ─> signal/<sample>_{plus,minus}.bw
              └─> genome BAM ─> RSeQC (infer_experiment, read_distribution,
                                       geneBody_coverage, inner_distance)

MultiQC <── FastQC (raw + trimmed), trimming reports, STAR Log.final.out,
            RSEM .stat/, RSeQC infer_experiment + read_distribution
  └─> qc/multiqc/multiqc_report.html
```

The Kallisto abundances, the gene body coverage and the inner distance files are
published but are not part of the MultiQC report; read those files directly.

### Stage Summary

| Stage | Tool | Input | Output | Reference |
|-------|------|-------|--------|-----------|
| 1. QC & Trimming | FastQC, Trim Galore | Raw FASTQ | Trimmed FASTQ + FastQC reports | references/01-qc-trimming.md |
| 2. Alignment | STAR (2-pass) | Trimmed FASTQ | Genome BAM + Transcriptome BAM + bedGraph | references/02-star-alignment.md |
| 3. Quantification | RSEM, Kallisto | Transcriptome BAM / trimmed FASTQ | Gene/transcript counts, TPM, FPKM | references/03-quantification.md |
| 4. Signal Tracks | bedGraphToBigWig | STAR bedGraph | Strand-specific bigWig | references/04-signal-tracks.md |
| 5. QC Metrics | RSeQC, MultiQC | Genome BAM, logs | Strandedness, read distribution, gene body coverage | references/05-qc-metrics.md |

## Input Requirements

### Required Files

- **RNA-seq FASTQ**: paired-end reads matched by the `--reads` glob (ENCODE standard;
  single-end with `--single_end`)
- **STAR genome index directory** (`--star_index`)
- **RSEM reference prefix** (`--rsem_index`) produced by `rsem-prepare-reference`
- **BED12 gene model for RSeQC** (`--rseqc_bed`)
- **Kallisto index file** (`--kallisto_index`) built with kallisto 0.50.1, unless
  `--skip_kallisto` is set. kallisto 0.50.1 writes index version 13 and rejects an index
  built with 0.48 or earlier

There is no sample sheet: samples are the pairs that `--reads` matches, and the sample
ID is the shared prefix of each pair. The gene annotation is not a workflow parameter
and there is no `--gtf`; the GTF is consumed when the STAR and RSEM references are built
(references/02 and references/03), so the annotation is fixed by the index you pass in.

## Library Strandedness

`--strandedness` takes one value for the whole run — `reverse` (default), `forward` or
`none` — and is validated before the first task. It drives three things at once: the RSEM
`--strandedness` flag, the kallisto strand flag, and which STAR bedGraph becomes which
signal track.

| Protocol | `--strandedness` | RSEM receives | Kallisto receives | Signal tracks |
|----------|------------------|---------------|-------------------|---------------|
| dUTP (ENCODE standard), Illumina TruSeq Stranded | `reverse` | `--strandedness reverse` | `--rf-stranded` | `<sample>_plus.bw`, `<sample>_minus.bw` |
| Directional ligation (some legacy protocols) | `forward` | `--strandedness forward` | `--fr-stranded` | `<sample>_plus.bw`, `<sample>_minus.bw` |
| SMARTer / SMART-Seq2 and other unstranded kits | `none` | `--strandedness none` | no strand flag | `<sample>_unstranded.bw` |

There is no per-sample strandedness and the workflow does not detect it. RSeQC
`infer_experiment.py` runs as a post-hoc check and writes
`qc/rseqc/<sample>.infer_experiment.txt`. If the library type is unknown, run a first
pass, read that file (references/05 explains the output), and rerun with the correct
`--strandedness` — the RSEM counts, the kallisto abundances and the signal tracks all
depend on it, so a wrong value has to be corrected by rerunning, not by post-processing.

## QC Thresholds

| Metric | Threshold | Produced by |
|--------|-----------|-------------|
| Total sequenced reads | >=30M PE reads | `fastqc/`, `star/<sample>.Log.final.out` |
| Uniquely mapped reads | >=70% of input reads | `star/<sample>.Log.final.out` |
| Multi-mapped reads | <10% | `star/<sample>.Log.final.out` |
| Strandedness agreement | >90% for a stranded library | `qc/rseqc/<sample>.infer_experiment.txt` |
| Exonic rate | >60% | `qc/rseqc/<sample>.read_distribution.txt` |
| Gene body coverage | Relatively uniform (5'/3' bias <1.5) | `qc/rseqc/<sample>.geneBody_coverage.geneBodyCoverage.txt` |

Not computed by this workflow: rRNA rate, library duplication rate (beyond the
sequence-level estimate inside the FastQC report), detected-gene counts, and saturation
curves. references/05-qc-metrics.md gives the commands to run those by hand on the
published BAM and RSEM output.

### Read Depth Guidelines

| Application | Minimum Reads (PE) | Recommended | Notes |
|-------------|-------------------|-------------|-------|
| Gene-level expression | 20M | 30M | ENCODE minimum |
| Transcript-level expression | 40M | 60M | Isoform resolution requires more depth |
| Differential expression | 20M per sample | 30M per sample | 3+ biological replicates per condition |
| Novel junction discovery | 60M | 100M+ | STAR 2-pass mode benefits from depth |
| Fusion detection | 50M | 80M+ | Chimeric reads are rare; needs a separate STAR run (references/02) |

## Execution

The versions the workflow runs are the ones in `scripts/Dockerfile`: STAR 2.7.11b,
RSEM 1.3.3, kallisto 0.50.1, samtools 1.19, RSeQC 5.0.3, Trim Galore 0.6.10, cutadapt 4.6,
MultiQC 1.21 and FastQC 0.12.1. The conda environment in `bioinformatics-installer`
(`environments/rnaseq-env.yml`) is a separate manual route pinned to the same versions of
STAR, RSEM, kallisto, samtools, RSeQC, Trim Galore, FastQC and MultiQC; it leaves cutadapt
to the Trim Galore package, adds salmon and subread, and does not carry
`bedGraphToBigWig`, which the signal-track step needs.

Every index flag is shown in the examples below because the defaults are bare names
resolved in the launch directory (`GRCh38_star_index`, `GRCh38_rsem_index/GRCh38`,
`gencode.v38.kallisto.idx`, `hg38_RefSeq.bed` for `--genome GRCh38`). The run stops
before the first task if any of them is missing. `--kallisto_index` is the one exception:
it is not read when `--skip_kallisto` is set.

### Quick Start (local, Docker)

```bash
nextflow run scripts/main.nf -profile local \
    --reads 'fastq/*_R{1,2}.fq.gz' \
    --genome GRCh38 \
    --star_index /ref/GRCh38_star_index \
    --rsem_index /ref/GRCh38_rsem_index/GRCh38 \
    --kallisto_index /ref/gencode.v38.kallisto.idx \
    --rseqc_bed /ref/hg38_RefSeq.bed \
    --strandedness reverse \
    --outdir results/
```

### SLURM HPC

```bash
nextflow run scripts/main.nf -profile slurm \
    --container /path/to/pipeline-rnaseq.sif \
    --slurm_queue normal \
    --reads 'fastq/*_R{1,2}.fq.gz' \
    --genome GRCh38 \
    --star_index /ref/GRCh38_star_index \
    --rsem_index /ref/GRCh38_rsem_index/GRCh38 \
    --kallisto_index /ref/gencode.v38.kallisto.idx \
    --rseqc_bed /ref/hg38_RefSeq.bed \
    --outdir results/
```

### Cloud

```bash
# Google Cloud Batch
nextflow run scripts/main.nf -profile gcp \
    --container us-docker.pkg.dev/<project>/<repo>/pipeline-rnaseq:1.0.0 \
    --gcp_project <project> \
    --gcp_workdir gs://<bucket>/work \
    --reads 'gs://<bucket>/fastq/*_R{1,2}.fq.gz' \
    --genome GRCh38 \
    --star_index gs://<bucket>/ref/GRCh38_star_index \
    --rsem_index gs://<bucket>/ref/GRCh38_rsem_index/GRCh38 \
    --kallisto_index gs://<bucket>/ref/gencode.v38.kallisto.idx \
    --rseqc_bed gs://<bucket>/ref/hg38_RefSeq.bed \
    --outdir gs://<bucket>/results

# AWS Batch
nextflow run scripts/main.nf -profile aws \
    --container <account>.dkr.ecr.<region>.amazonaws.com/pipeline-rnaseq:1.0.0 \
    --aws_queue <job-queue> \
    --aws_workdir s3://<bucket>/work \
    --reads 's3://<bucket>/fastq/*_R{1,2}.fq.gz' \
    --genome GRCh38 \
    --star_index s3://<bucket>/ref/GRCh38_star_index \
    --rsem_index s3://<bucket>/ref/GRCh38_rsem_index/GRCh38 \
    --kallisto_index s3://<bucket>/ref/gencode.v38.kallisto.idx \
    --rseqc_bed s3://<bucket>/ref/hg38_RefSeq.bed \
    --outdir s3://<bucket>/results
```

`--outdir` only sets where results are published; Google Batch and AWS Batch stage every
task through the work directory, and the workflow stops with an error if it or the
project/queue is missing.

## Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | required | Glob matching the FASTQ pairs, e.g. `'fastq/*_R{1,2}.fq.gz'`. Quote it |
| `--genome` | `GRCh38` | `GRCh38` or `mm10`; selects the default index names and nothing else |
| `--outdir` | `./results` | Directory results are published to |
| `--single_end` | `false` | Treat `--reads` as single files; kallisto then runs with the fixed `--single -l 200 -s 20`, and RSeQC `inner_distance.py` is skipped |
| `--strandedness` | `reverse` | `reverse`, `forward` or `none`; one value for the whole run |
| `--skip_kallisto` | `false` | Skip `KALLISTO_QUANT`; `--kallisto_index` is then not read |
| `--star_index` | `GRCh38_star_index` (`mm10_star_index`) | STAR genome directory |
| `--rsem_index` | `GRCh38_rsem_index/GRCh38` (`mm10_rsem_index/mm10`) | RSEM reference **prefix**, not a directory; every file starting with it is staged |
| `--kallisto_index` | `gencode.v38.kallisto.idx` (`gencode.vM27.kallisto.idx`) | kallisto index file; must be built with kallisto 0.50.1 (index version 13), not with 0.48 or earlier |
| `--rseqc_bed` | `hg38_RefSeq.bed` (`mm10_RefSeq.bed`) | BED12 gene model used by all four RSeQC modules |
| `--chrom_sizes` | `<star_index>/chrNameLength.txt` | Chromosome sizes for `bedGraphToBigWig` |

### Infrastructure parameters (`nextflow.config`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--container` | `encode-toolkit/pipeline-rnaseq:1.0.0` | Image built from `scripts/Dockerfile`. Pass a registry image for `gcp`/`aws`, or a `.sif` file for `slurm` |
| `--max_cpus`, `--max_memory`, `--max_time` | `16`, `64.GB`, `24.h` | Upper bounds applied to every process |
| `--slurm_queue`, `--slurm_account` | `normal`, none | SLURM partition and account |
| `--gcp_project`, `--gcp_workdir` | none (both required for `-profile gcp`) | Google Cloud project and `gs://` work directory |
| `--gcp_location`, `--gcp_disk` | `us-central1`, `200.GB` | Google Batch region and per-task disk |
| `--aws_queue`, `--aws_workdir` | none (both required for `-profile aws`) | AWS Batch job queue and `s3://` work directory |
| `--aws_region`, `--aws_cli_path` | `us-east-1`, `/home/ec2-user/miniconda/bin/aws` | AWS region, and the AWS CLI path inside the Batch AMI |

## Cloud Cost Estimates

| Platform | Instance | Cost/Sample | Time/Sample | Notes |
|----------|----------|-------------|-------------|-------|
| GCP | n1-highmem-8 | ~$3-6 | 2-4 hours | STAR index loading dominates; the `gcp` profile uses spot VMs |
| AWS | r5.2xlarge | ~$3-6 | 2-4 hours | r-series for STAR memory; spot recommended |
| Local | 8 cores, 36 GB | $0 | 3-6 hours | Docker required |
| SLURM | 8 cores, 36 GB | Varies | 2-4 hours | Singularity; pass the `.sif` with `--container` |

**Memory note**: `nextflow.config` asks for 36 GB for `STAR_ALIGN` (multiplied by the
attempt number on a retry, capped by `--max_memory`). The local executor refuses the task
on a machine with less, so a 32 GB host is not enough. If 36 GB is out of reach, rebuild the STAR index
with a larger `--genomeSAsparseD` (which shrinks the loaded index at some cost in
mapping speed) or move to a bigger machine. `--limitGenomeGenerateRAM` is a
`genomeGenerate` option and is not a parameter of this workflow.

## Output Directory Structure

```
results/
  fastqc/                                        # FastQC HTML/zip for raw and trimmed reads
  trimmed/
    <sample>_R1_val_1.fq.gz, <sample>_R2_val_2.fq.gz
    <sample>_R1.fq.gz_trimming_report.txt        # one per input file
  star/
    <sample>.Aligned.sortedByCoord.out.bam
    <sample>.Aligned.sortedByCoord.out.bam.bai
    <sample>.Aligned.toTranscriptome.out.bam     # RSEM input
    <sample>.Log.final.out
    <sample>.SJ.out.tab
    <sample>.ReadsPerGene.out.tab                # STAR gene counts
    <sample>.Signal.UniqueMultiple.str1.out.bg   # plus str2 unless --strandedness none
    <sample>.Signal.Unique.str1.out.bg           # same signal, unique mappers only
  rsem/
    <sample>.genes.results                       # gene_id, TPM, FPKM, expected_count
    <sample>.isoforms.results                    # transcript_id, TPM, FPKM, IsoPct
    <sample>.stat/                               # RSEM model statistics, read by MultiQC
  kallisto/                                      # absent with --skip_kallisto
    <sample>/abundance.tsv
    <sample>/abundance.h5                        # only from a kallisto build with HDF5 support
    <sample>/run_info.json
  signal/
    <sample>_plus.bw, <sample>_minus.bw          # or <sample>_unstranded.bw
  qc/
    rseqc/
      <sample>.infer_experiment.txt
      <sample>.read_distribution.txt
      <sample>.geneBody_coverage.*
      <sample>.inner_distance.*                  # paired-end only
    multiqc/
      multiqc_report.html
      multiqc_data/
  pipeline_info/
    timeline.html, report.html, trace.txt        # Nextflow execution reports
```

## Common Pitfalls

### 1. Insufficient Memory for STAR
STAR loads the whole genome index into memory. `STAR_ALIGN` is configured for 36 GB, so a
32 GB machine will not schedule the task at all under `-profile local`. On shared HPC
systems, check the per-job memory limit before submitting.

### 2. Wrong Strandedness Setting
Using incorrect strandedness results in near-zero gene counts. If you see uniformly low
counts, read `qc/rseqc/<sample>.infer_experiment.txt` and rerun with the matching
`--strandedness`. ENCODE dUTP libraries are `reverse` stranded.

### 3. Using FPKM for Cross-Sample Comparison
FPKM values are not comparable across samples because they depend on total library
composition. Use TPM (comparable across samples) or raw counts with DESeq2/edgeR
normalization for differential expression.

### 4. Ignoring Multi-Mapped Reads
RSEM uses an expectation-maximization algorithm to probabilistically assign multi-mapped
reads. This is critical for gene families and repetitive elements. Do not pre-filter
multi-mappers before RSEM quantification.

### 5. Assuming rRNA Contamination Was Checked
High rRNA contamination (>10%) indicates failed rRNA depletion and reduces effective
sequencing depth. This workflow does not measure it — run the manual check in
references/05-qc-metrics.md against `star/<sample>.Aligned.sortedByCoord.out.bam` before
trusting the quantifications.

### 6. Not Using 2-Pass Mode for Novel Junctions
STAR 1-pass mode only uses annotated splice junctions. 2-pass mode first discovers novel
junctions then re-maps, critical for non-model organisms or samples with extensive
alternative splicing. This workflow always runs `--twopassMode Basic`.

## Pipeline Scripts

| File | Description |
|------|-------------|
| `scripts/main.nf` | Nextflow DSL2 pipeline |
| `scripts/nextflow.config` | Execution profiles (local/slurm/gcp/aws) |
| `scripts/Dockerfile` | Docker build with STAR, RSEM, Kallisto, RSeQC |

## ENCODE Data Integration

After running on your own data, compare with ENCODE reference:

```python
# Find matching ENCODE RNA-seq experiments
encode_search_experiments(
    assay_title="total RNA-seq",
    organ="pancreas",
    biosample_type="tissue"
)

# Download ENCODE gene quantifications for comparison
encode_batch_download(
    download_dir="/data/encode_reference/",
    output_type="gene quantifications",
    assay_title="total RNA-seq",
    organ="pancreas",
    assembly="GRCh38"
)

# Download ENCODE signal tracks for browser visualization
encode_search_files(
    file_format="bigWig",
    assay_title="total RNA-seq",
    organ="pancreas",
    output_type="signal of unique reads"
)
```

## Pitfalls & Edge Cases

- **Strandedness must match library prep**: `--strandedness` feeds RSEM, kallisto and the
  signal tracks. Using the wrong value can halve gene counts, assign reads to antisense
  genes, and swap the plus/minus bigWigs.
- **rRNA contamination**: rRNA >10% wastes sequencing depth. Ribosomal depletion libraries
  should have <5%, poly-A selection libraries <1%. The workflow does not measure it; see
  references/05-qc-metrics.md for a manual count, or run Picard CollectRnaSeqMetrics
  outside the container.
- **STAR 2-pass mode is required**: The first pass discovers novel splice junctions; the
  second pass uses them. Single-pass STAR misses tissue-specific or rare splicing events,
  reducing sensitivity for differential exon usage.
- **Gene-level vs transcript-level quantification**: RSEM provides transcript-level estimates but gene-level aggregation is more robust for differential expression. Transcript-level analysis requires many more replicates (≥6).
- **TPM normalization is not for cross-sample comparison**: TPM normalizes within a sample but is NOT appropriate for comparing expression across conditions. Use DESeq2 size factors or TMM normalization for differential expression.
- **Batch effects in multi-lab data**: RNA-seq is highly sensitive to library prep method, sequencer, and lab. Always check for batch effects with PCA before combining datasets from different sources.

## Walkthrough: Processing ENCODE RNA-seq from FASTQ to Gene Quantification

**Goal**: Process raw RNA-seq FASTQ files through the ENCODE pipeline to generate gene expression quantifications (TPM/FPKM) and signal tracks.
**Context**: The ENCODE RNA-seq pipeline uses STAR 2-pass alignment and RSEM quantification, producing both gene-level and transcript-level expression estimates.

### Step 1: Find RNA-seq experiment

```
encode_get_experiment(accession="ENCSR000CPR")
```

Expected output:
```json
{
  "accession": "ENCSR000CPR",
  "assay_title": "total RNA-seq",
  "biosample_summary": "K562",
  "assembly": ["GRCh38"],
  "bio_replicate_count": 2,
  "status": "released"
}
```

### Step 2: List FASTQ files

```
encode_list_files(experiment_accession="ENCSR000CPR", file_format="fastq")
```

Expected output (a JSON array of file records; fields abridged):
```json
[
  {"accession": "ENCFF200RN1", "file_format": "fastq", "output_type": "reads", "file_size_human": "3.2 GB", "biological_replicates": [1], "status": "released"},
  {"accession": "ENCFF201RN2", "file_format": "fastq", "output_type": "reads", "file_size_human": "3.3 GB", "biological_replicates": [1], "status": "released"}
]
```

### Step 3: Download and name the FASTQs so a read-pair glob can find them

```
encode_download_files(file_accessions=["ENCFF200RN1", "ENCFF201RN2"], download_dir="/data/rnaseq/fastq")
```

ENCODE names every FASTQ after its accession (`ENCFF200RN1.fastq.gz`), with no `_R1`/`_R2`
in the name, so the two files of a pair share no prefix and the `--reads` glob cannot pair
them. Link them into the shape the glob expects. Which mate an accession is comes from the
ENCODE file record on encodeproject.org, which carries `paired_end` (1 or 2) and
`paired_with`; the MCP file tools do not return those two fields:

```bash
cd /data/rnaseq/fastq
ln -s ENCFF200RN1.fastq.gz k562_rep1_R1.fq.gz
ln -s ENCFF201RN2.fastq.gz k562_rep1_R2.fq.gz
```

### Step 4: Run the RNA-seq pipeline

```bash
nextflow run scripts/main.nf -profile local \
    --reads '/data/rnaseq/fastq/k562_*_R{1,2}.fq.gz' \
    --genome GRCh38 \
    --star_index /ref/GRCh38_star_index \
    --rsem_index /ref/GRCh38_rsem_index/GRCh38 \
    --kallisto_index /ref/gencode.v38.kallisto.idx \
    --rseqc_bed /ref/hg38_RefSeq.bed \
    --strandedness reverse \
    --outdir results/
```

Key pipeline steps:
1. FastQC on the raw reads, quality and adapter trimming with Trim Galore (which also runs FastQC on the trimmed reads)
2. STAR 2-pass alignment (splice-aware), writing the genome BAM, the transcriptome BAM, gene counts and bedGraphs
3. RSEM gene and isoform quantification (TPM, FPKM, expected counts)
4. Kallisto transcript quantification (optional, skipped with `--skip_kallisto`)
5. Signal track generation (bedGraph to bigWig)
6. RSeQC (infer_experiment, read_distribution, geneBody_coverage, inner_distance) and MultiQC

### Step 5: Validate output quality

| Metric | Threshold | Where to read it |
|---|---|---|
| Uniquely mapped rate | >= 70% | `star/<sample>.Log.final.out` |
| Strandedness agreement | > 90% and matching `--strandedness` | `qc/rseqc/<sample>.infer_experiment.txt` |
| Exonic rate | > 60% | `qc/rseqc/<sample>.read_distribution.txt` |
| Replicate correlation | >= 0.9 | compute yourself from the TPM column of `rsem/<sample>.genes.results` |

### Step 6: Use expression data with ENCODE epigenomic data

Compare gene expression with enhancer marks:
```
encode_search_experiments(assay_title="Histone ChIP-seq", biosample_term_name="K562", target="H3K27ac", organism="Homo sapiens")
```

**Interpretation**: Genes with high TPM AND nearby H3K27ac peaks have validated enhancer-gene connections. Low expression despite nearby enhancer marks suggests poised or tissue-specific regulation.

### Integration with downstream skills
- Gene quantifications feed into -> **peak-annotation** for expression-validated peak targets
- Expression data connects to -> **gtex-expression** for tissue comparison
- Processed data feeds into -> **compare-biosamples** for differential expression analysis
- Pipeline provenance logged by -> **data-provenance**

## Code Examples

### 1. Find RNA-seq experiments for a tissue

```
encode_search_experiments(assay_title="total RNA-seq", organ="liver", organism="Homo sapiens")
```

Expected output:
```json
{
  "results": [
    {"accession": "ENCSR300RNA", "assay_title": "total RNA-seq", "biosample_summary": "liver tissue male adult (54 years)", "status": "released"}
  ],
  "total": 35,
  "limit": 25,
  "offset": 0,
  "has_more": true,
  "next_offset": 25
}
```

### 2. Check for existing gene quantifications

```
encode_list_files(experiment_accession="ENCSR300RNA", file_format="tsv", output_type="gene quantifications", assembly="GRCh38")
```

Expected output (a JSON array of file records; fields abridged):
```json
[
  {"accession": "ENCFF400GEQ", "file_format": "tsv", "output_type": "gene quantifications", "assembly": "GRCh38", "file_size_human": "5.2 MB", "status": "released"}
]
```

### 3. Download expression data

```
encode_download_files(file_accessions=["ENCFF400GEQ"], download_dir="/data/rnaseq/quantification")
```

Expected output (fields abridged):
```json
{
  "downloaded": [
    {"accession": "ENCFF400GEQ", "file_path": "/data/rnaseq/quantification/ENCFF400GEQ.tsv", "file_size_human": "5.2 MB", "success": true, "md5_verified": true}
  ],
  "errors": [],
  "summary": {"total_requested": 1, "successful": 1, "failed": 0, "total_size_human": "5.2 MB"}
}
```

## Integration

| This skill produces... | Feed into... | Purpose |
|---|---|---|
| Gene expression (TPM/FPKM) | **peak-annotation** | Validate enhancer targets with expression data |
| Expression matrix | **gtex-expression** | Compare cell-line vs. tissue expression |
| Differential expression results | **compare-biosamples** | Identify tissue-specific gene regulation |
| Signal tracks (bigWig) | **visualization-workflow** | Display expression signal in genome browser |
| Expression quantifications | **disease-research** | Connect gene expression to disease phenotypes |
| Pipeline run parameters | **data-provenance** | Record STAR/RSEM versions and settings |
| QC metrics | **quality-assessment** | Validate against ENCODE RNA-seq standards |

## Related Skills

- **pipeline-guide** (parent): General pipeline selection and resource assessment
- **quality-assessment**: Deep-dive QC analysis beyond basic metrics
- **integrative-analysis**: Combine RNA-seq with ChIP-seq/ATAC-seq for regulatory inference
- **compare-biosamples**: Compare expression profiles across cell types
- **single-cell-encode**: For scRNA-seq data processing (different pipeline)
- **pipeline-chipseq**: Sibling pipeline for ChIP-seq data
- **pipeline-atacseq**: Sibling pipeline for ATAC-seq data
- **publication-trust**: Verify literature claims backing analytical decisions

## Presenting Results

When reporting RNA-seq pipeline results:

- **Mapping rate**: Report the STAR uniquely mapped rate (>70% expected), multi-mapped rate (<10%), and unmapped rate from `star/<sample>.Log.final.out`
- **Quantification paths**: Provide paths to `rsem/<sample>.genes.results` (TPM, FPKM, expected_count), `rsem/<sample>.isoforms.results`, and `kallisto/<sample>/abundance.tsv` when Kallisto ran
- **Strandedness**: Confirm the orientation reported in `qc/rseqc/<sample>.infer_experiment.txt` matches the `--strandedness` value the run used; if it does not, the run has to be repeated
- **Key QC metrics**: Present the exonic rate from `read_distribution.txt` and the gene body coverage uniformity from `geneBody_coverage.geneBodyCoverage.txt` in a summary table, alongside the MultiQC report at `qc/multiqc/multiqc_report.html`
- **Derived metrics**: Detected-gene counts (TPM>1), rRNA rate, library duplication rate and saturation are not produced by this workflow. Compute them separately if they are needed, and say so when reporting
- **Signal tracks**: Provide paths to `signal/<sample>_plus.bw` and `signal/<sample>_minus.bw` (or `signal/<sample>_unstranded.bw` for an unstranded run)
- **Next steps**: Suggest `integrative-analysis` to combine RNA-seq with ChIP-seq/ATAC-seq for regulatory inference, or `compare-biosamples` for cross-tissue expression comparison

## For the request: "$ARGUMENTS"
