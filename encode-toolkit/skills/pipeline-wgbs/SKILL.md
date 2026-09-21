---
name: pipeline-wgbs
description: "Execute ENCODE Whole Genome Bisulfite Sequencing (WGBS) pipeline from FASTQ to methylation calls. Child of pipeline-guide. Provides Nextflow execution with Docker and cloud deployment. Use when processing WGBS/bisulfite-seq data, calling methylation levels, generating bedMethyl files. Trigger on: WGBS pipeline, bisulfite sequencing, methylation calling, DNA methylation pipeline, bismark, bwa-meth, bedMethyl."
---

# ENCODE WGBS Pipeline: FASTQ to Methylation Calls

## When to Use

- User wants to run a WGBS/bisulfite sequencing pipeline from FASTQ to methylation calls
- User asks about "WGBS pipeline", "bisulfite sequencing", "methylation calling", "Bismark", or "bedMethyl"
- User needs to process whole-genome bisulfite sequencing data following ENCODE standards
- Example queries: "process my WGBS FASTQs", "call methylation levels from bisulfite-seq", "run Bismark on my WGBS data"

Execute the ENCODE DNA methylation pipeline for Whole Genome Bisulfite Sequencing data,
producing per-CpG methylation levels in bedMethyl format.

## Pipeline Overview

```
FASTQ -> Trim adapters -> Bismark align -> Deduplicate -> sort + index -> MethylDackel -> bedMethyl
  |           |                 |               |                |             |             |
 FastQC   Trim Galore    Bismark (Bowtie2)  deduplicate_   samtools      mbias +       coverage
                                             bismark                     extract        stats
```

Paired-end only: `TRIM_GALORE` and `BISMARK_ALIGN` both take two mates, and `--reads` must
match pairs.

### ENCODE Repository

- **GitHub**: `ENCODE-DCC/dna-me-pipeline`
- **Container**: built from `scripts/Dockerfile` in this skill (`docker build -t encode-toolkit/pipeline-wgbs:1.0.0 scripts/`); override with `--container`
- **WDL**: Available for Cromwell execution
- **This skill**: Nextflow DSL2 reimplementation for portability

## Core Tools and Versions

These are the versions in `scripts/Dockerfile`, which is what the workflow runs.

| Tool | Version | Purpose | Citation |
|------|---------|---------|----------|
| Trim Galore | 0.6.10 | Adapter + quality trimming (bisulfite-aware) | Krueger (Babraham) |
| Bismark | 0.24.2 | Bisulfite-aware alignment + deduplication | Krueger & Andrews 2011 |
| Bowtie2 | 2.5.4 | Backend aligner used by Bismark | Langmead & Salzberg 2012 |
| MethylDackel | 0.6.1 | Methylation extraction from BAM | Ryan (GitHub) |
| samtools | 1.19 | BAM sorting and indexing | Li et al. 2009 |
| htslib | 1.19 | `bgzip` and `tabix` for the bedMethyl files | Bonfield et al. 2021 |
| FastQC | 0.12.1 | Read quality assessment | Andrews (Babraham) |
| MultiQC | 1.21 | Aggregated QC reporting | Ewels et al. 2016 |

The conda environment in `bioinformatics-installer` (`environments/wgbs-env.yml`) is a
separate manual route pinned to the same versions of these tools; it also carries
bedtools, which is not in the image.

## Key Literature

1. **Krueger & Andrews 2011** - "Bismark: a flexible aligner and methylation caller for
   Bisulfite-Seq applications" (Bioinformatics, ~4,000 citations)
   DOI: 10.1093/bioinformatics/btr167

2. **Lister et al. 2009** - "Human DNA methylomes at base resolution show widespread
   epigenomic differences" (Nature, ~5,000 citations)
   DOI: 10.1038/nature08514

3. **Schultz et al. 2015** - "Human body epigenome maps reveal noncanonical DNA
   methylation variation" (Nature, ~1,500 citations)
   DOI: 10.1038/nature14248

4. **Pedersen et al. 2014** - "Fast and accurate alignment of long bisulfite-seq reads"
   arXiv:1401.1129 (bwa-meth)

5. **Amemiya et al. 2019** - "The ENCODE Blacklist" (Scientific Reports, ~1,372 citations)
   DOI: 10.1038/s41598-019-45839-z

## Execution

`--reads` and `--genome_dir` are both required; the workflow stops before the first task
if either is missing. `--genome_dir` is a Bismark genome folder — the output of
`bismark_genome_preparation`, which must still contain the genome `.fa`, because
MethylDackel reads it from there.

### Quick Start (local, Docker)

```bash
nextflow run scripts/main.nf -profile local \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --genome_dir /ref/bismark_index \
    --min_coverage 5 \
    --outdir results/ \
    -resume
```

### SLURM HPC

```bash
nextflow run scripts/main.nf -profile slurm \
    --container /path/to/pipeline-wgbs.sif \
    --slurm_queue normal \
    --reads '/data/fastq/*_R{1,2}.fastq.gz' \
    --genome_dir /ref/bismark_index \
    --outdir results/ \
    -resume
```

### Cloud

```bash
# Google Cloud Batch
nextflow run scripts/main.nf -profile gcp \
    --container us-docker.pkg.dev/<project>/<repo>/pipeline-wgbs:1.0.0 \
    --gcp_project <project> \
    --gcp_workdir gs://<bucket>/work \
    --reads 'gs://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --genome_dir gs://<bucket>/ref/bismark_index \
    --outdir gs://<bucket>/results

# AWS Batch
nextflow run scripts/main.nf -profile aws \
    --container <account>.dkr.ecr.<region>.amazonaws.com/pipeline-wgbs:1.0.0 \
    --aws_queue <job-queue> \
    --aws_workdir s3://<bucket>/work \
    --reads 's3://<bucket>/fastq/*_R{1,2}.fastq.gz' \
    --genome_dir s3://<bucket>/ref/bismark_index \
    --outdir s3://<bucket>/results
```

`--outdir` only sets where results are published; Google Batch and AWS Batch stage every
task through the work directory, and the workflow stops with an error if it or the
project/queue is missing.

## Resource Requirements

| Step | CPUs | RAM | Time (30x human) |
|------|------|-----|-------------------|
| FastQC | 2 | 4 GB | <1 hour |
| Trim Galore | 4 | 4 GB | 1-2 hours |
| Bismark align | 8 | 48 GB | 8-16 hours |
| Deduplication | 2 | 16 GB | 1-2 hours |
| Sort + index | 4 | 8 GB | 1-2 hours |
| MethylDackel (mbias, extract) | 2, 4 | 8 GB | 1-2 hours |
| **Total** | **8** | **48 GB** | **12-24 hours** |

Every process asks for `memory { N.GB * task.attempt }`, so a task killed for running out
of memory is retried with more: the second attempt gets twice the figure in the table, the
third three times it, bounded by `--max_memory`. `nextflow.config` scales the time the
same way for `BISMARK_ALIGN`, `DEDUPLICATE` and `METHYLDACKEL_EXTRACT`, bounded by
`--max_time`; the other processes declare no time limit. A task is retried only for exit
codes 130-145 and 104 (killed for exceeding a limit); any other failure stops the run.

## Pipeline Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | required | Glob matching the paired FASTQ files, e.g. `'/data/fastq/*_R{1,2}.fastq.gz'`. Quote it |
| `--genome_dir` | required | Bismark genome folder (`bismark_genome_preparation` output, including the genome `.fa`) |
| `--outdir` | `./results` | Directory results are published to |
| `--min_coverage` | `5` | Minimum read count a site must reach to appear in the bedMethyl files. The MethylDackel bedGraphs are written unfiltered |
| `--merge_context` | `true` | Merge the two strands of each CpG/CHG into one record (`MethylDackel --mergeContext`) |
| `--skip_dedup` | `false` | Skip `deduplicate_bismark`; the sorted BAM is then the raw alignment |

### Infrastructure parameters (`nextflow.config`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--container` | `encode-toolkit/pipeline-wgbs:1.0.0` | Image built from `scripts/Dockerfile`. Pass a registry image for `gcp`/`aws`, or a `.sif` file for `slurm` |
| `--max_cpus`, `--max_memory`, `--max_time` | `16`, `64.GB`, `48.h` | Upper bounds applied to every process |
| `--slurm_queue`, `--slurm_account` | `normal`, none | SLURM partition and account |
| `--gcp_project`, `--gcp_workdir` | none (both required for `-profile gcp`) | Google Cloud project and `gs://` work directory |
| `--gcp_location`, `--gcp_disk` | `us-central1`, `200.GB` | Google Batch region and per-task disk |
| `--aws_queue`, `--aws_workdir` | none (both required for `-profile aws`) | AWS Batch job queue and `s3://` work directory |
| `--aws_region`, `--aws_cli_path` | `us-east-1`, `/home/ec2-user/miniconda/bin/aws` | AWS region, and the AWS CLI path inside the Batch AMI |

Notes on what the workflow does and does not do:
- **Aligner**: Bismark with Bowtie2. bwa-meth is described in `references/02-bismark-alignment.md`
  as a manual alternative; it is not an option of this workflow and is not in the image.
- **Deduplication**: `deduplicate_bismark`. Picard is not installed and is not used.
- **BAM filtering**: there is no MAPQ or flag filtering step. The BAM handed to
  MethylDackel is the Bismark output, deduplicated and sorted. `references/03-dedup-filtering.md`
  gives the manual filtering commands if you want them.
- **Overlapping mates**: MethylDackel never counts both mates of an overlapping read pair, so
  there is no switch for it. `--merge_context` is a separate choice about per-CpG versus
  per-cytosine output.
- **Read-position trimming at extraction**: none. `MethylDackel extract` is called without
  `--nOT`/`--nOB`, because `--clip_R2 10` at the trimming step already removes the
  end-repair bias. See the M-bias pitfall below.
- **Bisulfite conversion rate**: not computed by the workflow. Estimate it from the CHH
  methylation percentage in the Bismark alignment report, or align a lambda/pUC19 spike-in
  separately when the sample carries real non-CpG methylation (ESCs, neurons).
  `references/05-qc-metrics.md` has both commands.

## Output Files

```
results/
  fastqc/                      # FastQC on the raw reads
  trim_galore/                 # trimmed reads, trimming reports, and FastQC on the trimmed reads
  bismark/
    alignments/
      {sample}.sorted.bam      # final BAM: deduplicated unless --skip_dedup, sorted
      {sample}.sorted.bam.bai
      *_PE_report.txt          # Bismark alignment report (mapping rate, context methylation)
    dedup_reports/
      *.deduplication_report.txt
    mbias/
      {sample}_mbias_*.svg     # M-bias plots, one per strand
      {sample}_mbias_report.txt  # MethylDackel's suggested inclusion bounds
    methylation/
      {sample}_CpG.bedGraph    # MethylDackel output, every covered cytosine/CpG
      {sample}_CHG.bedGraph
      {sample}_CHH.bedGraph
      {sample}.CpG.bedMethyl.gz (+ .tbi)   # primary output, filtered by --min_coverage
      {sample}.CHG.bedMethyl.gz (+ .tbi)   # non-CpG contexts
      {sample}.CHH.bedMethyl.gz (+ .tbi)
  coverage/
    {sample}.coverage_stats.txt
  multiqc/
    multiqc_report.html
  pipeline_info/
    timeline.html, report.html, trace.txt  # Nextflow execution reports
```

Note the naming: the MethylDackel bedGraphs use an underscore before the context
(`{sample}_CpG.bedGraph`), the bedMethyl files a dot (`{sample}.CpG.bedMethyl.gz`).

### bedMethyl Format

The primary output is per-CpG methylation in the ENCODE bedMethyl layout
(https://www.encodeproject.org/data-standards/wgbs/):

```
chr1  10468  10470  .  12  .  10468  10470  0,0,0  12  83
```

Columns: chr, start, end, name (`.`), score, strand, thickStart, thickEnd, colour,
coverage, percent methylated. The score in column 5 is the read count capped at 1000, the
strand in column 6 is `.` because MethylDackel's bedGraph carries no strand, and column 11
is an integer percentage. Sites below `--min_coverage` are dropped.

### Coverage Statistics

`coverage/{sample}.coverage_stats.txt` is computed from the unfiltered
`{sample}_CpG.bedGraph`, so it describes every covered site, not just the reported ones.
It has five lines — this is a run at the default `--min_coverage 5`:

```
Covered CpGs (>=1x): 27184023
Mean coverage of covered CpGs: 12.4
Covered CpGs >=5x: 22903511 (84.3%)
Covered CpGs >=10x: 16992841 (62.5%)
Covered CpGs >=5x (--min_coverage, kept in bedMethyl): 22903511 (84.3%)
```

The 5x and 10x lines are fixed thresholds and do not follow `--min_coverage`. The last
line does, and it is the one that describes what actually reached the bedMethyl files: at
the default it repeats the >=5x line, with `--min_coverage 10` it would repeat the >=10x
line, and at any other value it stands on its own. All the percentages are of *covered*
CpGs, not of all CpGs in the genome — the workflow never counts genomic CpGs that got zero
reads. With `--merge_context false` the records are per cytosine rather than per CpG.

## QC Thresholds (ENCODE Standards)

| Metric | Pass | Warning | Fail | Source |
|--------|------|---------|------|--------|
| Mapping rate | >70% | 50-70% | <50% | `bismark/alignments/*_PE_report.txt` |
| Duplication rate | <30% | 30-50% | >50% | `bismark/dedup_reports/*.deduplication_report.txt` |
| Mean coverage of covered CpGs | >10x | 5-10x | <5x | `coverage/{sample}.coverage_stats.txt` |
| Covered CpGs reaching >=5x | >80% | 60-80% | <60% | `coverage/{sample}.coverage_stats.txt` |
| Bisulfite conversion rate | ≥98% | 95-98% | <95% | External; see `references/05-qc-metrics.md` |
| Lambda spike-in conversion | ≥98% | 95-98% | <95% | External; requires a separate alignment |

## Critical Pitfalls

### RRBS vs WGBS
RRBS (Reduced Representation) uses MspI digestion and covers ~10% of CpGs.
WGBS covers the full genome. These are DIFFERENT protocols:
- RRBS: this workflow is not set up for it. Deduplication can be skipped
  (`--skip_dedup true`), but the trimming is WGBS-specific and hardcoded
  (`--clip_R2 10 --three_prime_clip_R1 1`), there is no `--rrbs` switch, and the workflow
  re-trims whatever FASTQs it is given. Run RRBS trimming and alignment by hand —
  `references/01-qc-trimming.md` has the command
- WGBS: full dedup required, standard Trim Galore settings
- Never mix RRBS and WGBS data in the same analysis

### Strand-Specific vs Merged CpG
Bismark reports methylation per strand by default. For most analyses, merge
complementary CpG strands:
- Forward C at position N and reverse G at position N+1 are the same CpG
- MethylDackel `--mergeContext` handles this automatically
- `--merge_context true` is the workflow default; set it to `false` only when you need
  strand-specific data

### Incomplete Bisulfite Conversion
Conversion artifacts produce false methylation calls:
- Always include lambda phage or pUC19 spike-in DNA
- Unmethylated spike-in should show ≥98% conversion
- The workflow does not measure this. Align the spike-in separately, or read the CHH
  methylation percentage out of the Bismark report, before trusting the calls

### M-bias Plots
`MethylDackel mbias` runs on every sample and writes plots plus
`bismark/mbias/{sample}_mbias_report.txt`, which captures MethylDackel's suggested
inclusion bounds. Extraction itself excludes no read positions:
- End-repair artifacts cause elevated methylation at read ends. The 5' end-repair bias of
  read 2 is already removed at the trimming step by `--clip_R2 10`
- If the plots still show bias, re-run `MethylDackel extract` manually with `--OT`/`--OB`
  set to the suggested bounds. Note the argument order: `--OT A,B,C,D` *includes*
  positions A-B on read 1 and C-D on read 2, and the matching `--nOT a,b,c,d` *excludes*
  a bases from the start of read 1, b from its end, c from the start of read 2 and d from
  its end
- `--maxDepth` is not an option of MethylDackel 0.6.1 and must not be passed

### Low Coverage Regions
Regions with <5x coverage have unreliable methylation estimates:
- The bedMethyl files are filtered to `--min_coverage` (default 5); the bedGraphs are not
- For differential methylation analysis, consider `--min_coverage 10`
- Report the >=5x and >=10x fractions from `coverage/{sample}.coverage_stats.txt`. Those
  two thresholds are fixed in the workflow and do not follow `--min_coverage`; the last
  line of the file is the one that does, so quote it too whenever the run used a value
  other than 5

## Provenance Integration

After pipeline completion, log all outputs:

```python
# Log derived bedMethyl files
encode_log_derived_file(
    file_path="/results/bismark/methylation/sample1.CpG.bedMethyl.gz",
    source_accessions=["ENCSR...", "ENCFF..."],
    description="CpG methylation calls from ENCODE WGBS pipeline",
    file_type="bedMethyl",
    tool_used="Bismark 0.24.2 + MethylDackel 0.6.1",
    parameters="bismark --genome /ref -1 R1.fq.gz -2 R2.fq.gz; MethylDackel extract --mergeContext --CHG --CHH; bedMethyl filtered at --min_coverage 5"
)
```

## Reference Files

Detailed step-by-step documentation is provided in the `references/` directory:

1. `01-qc-trimming.md` -- Bisulfite-specific adapter trimming with Trim Galore
2. `02-bismark-alignment.md` -- Bismark alignment and the bwa-meth manual alternative
3. `03-dedup-filtering.md` -- Deduplication, and manual BAM filtering the workflow skips
4. `04-methylation-calling.md` -- MethylDackel extraction and bedMethyl generation
5. `05-qc-metrics.md` -- Coverage stats, M-bias, and manual conversion-rate checks

## Walkthrough: Processing ENCODE WGBS from FASTQ to Methylation Calls

**Goal**: Process whole-genome bisulfite sequencing FASTQ files through the ENCODE pipeline to generate per-CpG methylation calls for epigenomic analysis.
**Context**: WGBS requires bisulfite-aware alignment (Bismark) and per-CpG methylation extraction (MethylDackel), with ≥98% bisulfite conversion expected of the library.

### Step 1: Find WGBS experiment

```
encode_get_experiment(accession="ENCSR765JPC")
```

Expected output:
```json
{
  "accession": "ENCSR765JPC",
  "assay_title": "WGBS",
  "biosample_summary": "liver tissue male adult (54 years)",
  "assembly": ["GRCh38"],
  "bio_replicate_count": 2,
  "tech_replicate_count": 2,
  "status": "released"
}
```

### Step 2: List FASTQ files

```
encode_list_files(experiment_accession="ENCSR765JPC", file_format="fastq")
```

Expected output (a JSON array of file records; fields abridged):
```json
[
  {"accession": "ENCFF300BS1", "file_format": "fastq", "output_type": "reads", "file_size_human": "45.0 GB", "biological_replicates": [1], "status": "released"},
  {"accession": "ENCFF301BS2", "file_format": "fastq", "output_type": "reads", "file_size_human": "46.0 GB", "biological_replicates": [1], "status": "released"}
]
```

**Interpretation**: WGBS files are very large (~45GB per read file). Ensure adequate storage (>500GB for processing).

### Step 3: Download and name the FASTQs so a read-pair glob can find them

```
encode_download_files(file_accessions=["ENCFF300BS1", "ENCFF301BS2"], download_dir="/data/wgbs/fastq")
```

ENCODE names every FASTQ after its accession (`ENCFF300BS1.fastq.gz`), with no `_R1`/`_R2`
in the name, so the two files of a pair share no prefix and the `--reads` glob cannot pair
them. Link them into the shape the glob expects. Which mate an accession is comes from the
ENCODE file record on encodeproject.org, which carries `paired_end` (1 or 2) and
`paired_with`; the MCP file tools do not return those two fields:

```bash
cd /data/wgbs/fastq
ln -s ENCFF300BS1.fastq.gz liver_rep1_R1.fastq.gz
ln -s ENCFF301BS2.fastq.gz liver_rep1_R2.fastq.gz
```

### Step 4: Run the WGBS pipeline

```bash
nextflow run scripts/main.nf -profile local \
    --reads '/data/wgbs/fastq/liver_*_R{1,2}.fastq.gz' \
    --genome_dir /ref/bismark_index \
    --min_coverage 5 \
    --outdir results/ \
    -resume
```

Key pipeline steps:
1. FastQC on the raw reads, then adapter/quality trimming with Trim Galore (which also runs FastQC on the trimmed reads)
2. Bisulfite-aware alignment (Bismark with Bowtie2)
3. Deduplication (`deduplicate_bismark`), then sort and index
4. M-bias assessment and methylation extraction (MethylDackel), then conversion to bedMethyl
5. CpG coverage statistics
6. MultiQC aggregation

### Step 5: Validate output quality

| Metric | Threshold | Purpose |
|---|---|---|
| Mapping rate | > 70% | Alignment success (Bismark report) |
| Duplication rate | < 30% | Library complexity (dedup report) |
| Covered CpGs at >= 10x | Enough for DMR calling | Statistical power (coverage stats) |
| Bisulfite conversion | >= 98% | Library quality; measured outside this workflow |

### Step 6: Identify differentially methylated regions

Feed per-CpG methylation into -> **methylation-aggregation** for cross-tissue comparison and HMR/UMR/PMD identification.

### Integration with downstream skills
- Per-CpG methylation files feed into -> **methylation-aggregation** for cross-tissue atlas
- DMRs feed into -> **peak-annotation** for nearest gene assignment
- Methylation at regulatory elements connects to -> **regulatory-elements**
- CpG variant methylation integrates with -> **variant-annotation**
- Pipeline provenance logged by -> **data-provenance**

## Code Examples

### 1. Find WGBS data for methylation analysis

```
encode_search_experiments(
  assay_title="WGBS",
  organ="brain"
)
```

Expected output:
```json
{
  "results": [
    {
      "accession": "ENCSR321BRN",
      "assay_title": "WGBS",
      "biosample_summary": "brain tissue female adult (53 years)",
      "organ": "brain",
      "status": "released"
    }
  ],
  "total": 6,
  "limit": 25,
  "offset": 0,
  "has_more": false,
  "next_offset": null
}
```

### 2. Download processed methylation files

```
encode_search_files(
  assay_title="WGBS",
  organ="brain",
  file_format="bed",
  output_type="methylation state at CpG",
  assembly="GRCh38"
)
```

Expected output (fields abridged):
```json
{
  "results": [
    {
      "accession": "ENCFF567MET",
      "file_format": "bed",
      "output_type": "methylation state at CpG",
      "assembly": "GRCh38",
      "file_size": 886144860,
      "file_size_human": "845.1 MB",
      "experiment_accession": "ENCSR321BRN"
    }
  ],
  "total": 3,
  "limit": 25,
  "offset": 0,
  "has_more": false,
  "next_offset": null
}
```

## Integration

| This skill produces... | Feed into... | Purpose |
|---|---|---|
| Per-CpG methylation (bedMethyl) | **methylation-aggregation** | Cross-tissue methylation atlas |
| Per-CpG methylation (bedMethyl) | **visualization-workflow** | Display methylation in a genome browser (tabix-indexed BED; the workflow writes no bigWig) |
| Differentially methylated regions | **peak-annotation** | Assign DMRs to nearest genes |
| Methylation at regulatory sites | **regulatory-elements** | Correlate methylation with cCRE activity |
| CpG methylation near variants | **variant-annotation** | Annotate variants affecting CpG methylation |
| Coverage and mapping metrics | **quality-assessment** | Validate coverage and alignment against ENCODE standards |
| Pipeline parameters | **data-provenance** | Record Bismark/MethylDackel versions |
| Methylation at promoters | **gtex-expression** | Correlate promoter methylation with gene expression |

## Related Skills

- `pipeline-guide` -- Parent skill with compute resource assessment and cloud setup
- `methylation-aggregation` -- Aggregate methylation data across samples/tissues
- `quality-assessment` -- Evaluate pipeline output quality metrics
- `data-provenance` -- Track all pipeline inputs, outputs, and parameters
- `download-encode` -- Download ENCODE WGBS FASTQ files for pipeline input
- `publication-trust` -- Verify literature claims backing analytical decisions

## Presenting Results

When reporting WGBS pipeline results:

- **Mapping and duplication**: Report the mapping rate from `bismark/alignments/*_PE_report.txt` (>70% expected) and the duplication rate from `bismark/dedup_reports/` (<30% expected)
- **CpG coverage depth**: Report the mean coverage of covered CpGs, the percentage reaching >=5x and >=10x, and the last line's `--min_coverage` count from `coverage/{sample}.coverage_stats.txt`, and state that these are percentages of covered CpGs rather than of all genomic CpGs
- **Global methylation level**: Report the genome-wide average CpG methylation percentage and note any non-CpG (CHG/CHH) methylation if relevant to the tissue. Both come from the Bismark report or from the bedGraphs
- **bedMethyl output paths**: Provide paths to the CpG bedMethyl file and the CHG/CHH context files, and state the `--min_coverage` value they were filtered at
- **M-bias assessment**: No read positions are excluded at extraction. Report MethylDackel's suggested bounds from `bismark/mbias/{sample}_mbias_report.txt`, and say whether a manual re-extraction with `--OT`/`--OB` was done
- **Bisulfite conversion rate**: This workflow does not compute it. If a lambda/pUC19 spike-in was aligned separately, or the CHH percentage was used as a proxy, report the number and say which method produced it (≥98% pass, 95-98% warning, <95% fail)
- **Output summary**: Report covered CpG sites, sites passing the coverage filter, and the MultiQC report at `multiqc/multiqc_report.html`
- **Next steps**: Suggest `methylation-aggregation` for cross-sample averaging and HMR/UMR/PMD identification

## For the request: "$ARGUMENTS"
