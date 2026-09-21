---
name: data-provenance
description: Track exact provenance for every operation on ENCODE data — tool versions, reference files, scripts, parameters, and timestamps — to enable publication-ready methods writing. Use when the user processes ENCODE files, runs any bioinformatics tool, creates filtered/merged datasets, runs pipelines, performs liftover, uses R/Python/Bash for analysis, or needs to document their analysis chain for reproducibility and publication. Also use when the user says "write me methods" to auto-generate methods sections from the provenance log. This skill implements comprehensive provenance documentation: every tool, every version, every reference file, every parameter, every accession — no shortcuts. Use this skill for ANY processing step, ANY file transformation, ANY analysis operation on ENCODE data.
---

# Exact Provenance Tracking and Methods Writing

## When to Use

- User wants to track the full analysis chain from ENCODE download through processing to publication figure
- User asks about "provenance", "reproducibility", "methods section", or "analysis log"
- User needs to log derived files with their processing parameters for audit trail
- User wants to auto-generate publication-ready methods text from logged analysis steps
- Example queries: "log my peak calling step", "generate a methods section from my analysis", "show the provenance chain for this figure"

Track every operation on ENCODE data with exact tool versions, reference files, scripts, parameters, and timestamps to enable publication-ready methods sections.

## Scientific Rationale

**The question**: "What exactly was done to this data, and can someone else reproduce it identically?"

Reproducibility is the foundation of science. Yet the "Methods" sections of most genomics papers are vague — "reads were aligned with STAR" tells you nothing about which STAR version, which genome index, which parameters, or which annotation version was used. The difference between GENCODE v38 and v39 gene annotations can change thousands of gene assignments.

### Comprehensive Provenance Standard

This skill implements a documentation standard where every operation records:
1. **Tool**: Exact name and version (e.g., `bedtools v2.31.0`, not just "bedtools")
2. **Reference files**: Exact source, version, and download URL (e.g., "GRCh38.p14 chromosome sizes from UCSC, downloaded 2024-01-15")
3. **Parameters**: Complete command or function call, not just key parameters
4. **Input**: Exact file paths, accessions, and checksums
5. **Output**: File path, description, and checksum
6. **Script**: If a custom script was used, store the script alongside the output
7. **Timestamp**: When the operation was performed
8. **Environment**: R version, Python version, package versions, OS

This creates a complete audit trail such that a methods section can be auto-generated with zero ambiguity.

### Why This Level of Detail Matters

Consider a simple liftover operation. A vague log says "coordinates were lifted from hg19 to hg38." A comprehensive provenance log says:

> "Genomic coordinates were lifted from GRCh37/hg19 to GRCh38/hg38 using UCSC liftOver (v377, Kent et al. 2002, PMID: 12045153). The chain file hg19ToHg38.over.chain.gz was obtained from UCSC Genome Browser (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/, accessed 2024-01-15, MD5: abc123...). Of 45,231 input regions, 44,892 (99.25%) were successfully converted; 339 regions (0.75%) failed to map and were excluded. Unmapped regions were logged to unmapped.bed."

The second version can be reproduced exactly. The first cannot.

## Step 1: Initialize Experiment Log

At the start of any analysis session, create an experiment log:

### Log Structure
```
project_dir/
├── experiment_log.json          # Machine-readable provenance log
├── scripts/                     # All scripts used in this analysis
│   ├── 001_download.sh
│   ├── 002_filter_peaks.sh
│   └── 003_merge_samples.R
├── reference_files/             # Reference files used (or symlinks)
│   ├── GRCh38.chrom.sizes
│   └── gencode.v44.annotation.gtf
├── data/                        # ENCODE downloads
│   └── (organized by experiment)
├── derived/                     # All derived files
│   ├── filtered_peaks/
│   └── merged_results/
└── methods/                     # Auto-generated methods text
    └── methods_draft.md
```

### Experiment Log Format (experiment_log.json)
```json
{
  "project": "H3K27ac analysis in human pancreas",
  "created": "2024-01-15T10:30:00Z",
  "analyst": "Dr. A. Mawla",
  "organism": "Homo sapiens",
  "assembly": "GRCh38",
  "gene_annotation": "GENCODE v44",
  "operations": [],
  "encode_experiments": [],
  "software_environment": {},
  "reference_files": []
}
```

### Track Source ENCODE Experiments
```
encode_track_experiment(accession="ENCSR...", notes="Experiment log entry")
```

For each experiment, record in the log:
| Field | Example | Source |
|-------|---------|--------|
| Accession | ENCSR133RZO | ENCODE portal |
| Assay | Histone ChIP-seq | encode_get_experiment |
| Target | H3K27ac | encode_get_experiment |
| Biosample | pancreas tissue | encode_get_experiment |
| Lab | Bing Ren, UCSD | encode_get_experiment |
| Replicates | 2 biological | encode_get_experiment (`bio_replicate_count`) |
| Sequencer | Illumina HiSeq 4000 | ENCODE portal |
| Read length | 76bp PE | ENCODE portal |
| Read count | 42.3M per rep | File metadata |
| Library | TruSeq ChIP | ENCODE portal |
| Batch/date | 2019-06-15 | encode_get_experiment |

## Step 2: Log Every Operation

### Operation Log Entry Format
Every processing step creates a log entry with these fields:

```json
{
  "operation_id": "op_003",
  "timestamp": "2024-01-15T14:22:00Z",
  "description": "Filter H3K27ac peaks by signalValue",
  "category": "filtering",
  "tool": {
    "name": "bedtools",
    "version": "2.31.0",
    "citation": "Quinlan & Hall 2010, Bioinformatics, DOI:10.1093/bioinformatics/btq033"
  },
  "command": "awk '$7 >= 4.5' ENCFF123ABC.bed | bedtools intersect -a stdin -b blacklist.bed -v > filtered_peaks.bed",
  "script_path": "scripts/002_filter_peaks.sh",
  "inputs": [
    {
      "file": "ENCFF123ABC.bed",
      "accession": "ENCFF123ABC",
      "type": "IDR thresholded peaks",
      "md5": "abc123..."
    }
  ],
  "reference_files": [
    {
      "file": "hg38-blacklist.v2.bed.gz",
      "source": "ENCODE Blacklist v2 (Amemiya et al. 2019, Sci Rep, DOI:10.1038/s41598-019-45839-z)",
      "url": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
      "md5": "def456..."
    }
  ],
  "parameters": {
    "signalValue_threshold": 4.5,
    "blacklist_filter": "exclude overlapping regions"
  },
  "outputs": [
    {
      "file": "derived/filtered_peaks/H3K27ac_pancreas_filtered.bed",
      "description": "H3K27ac peaks in pancreas, signalValue >= 4.5, blacklist-filtered",
      "regions_count": 34521,
      "md5": "ghi789..."
    }
  ],
  "statistics": {
    "input_regions": 45231,
    "output_regions": 34521,
    "filtered_out": 10710,
    "filter_rate": "23.7%"
  }
}
```

### Common Operations to Log

#### Downloading ENCODE Files
```
encode_download_files(
    file_accessions=["ENCFF..."],
    download_dir="/path/to/data/",
    organize_by="experiment",
    verify_md5=True
)
```
Log: file accession, download URL, MD5 verification result, file size, download timestamp.

#### Genome Coordinate Liftover
Log: liftOver version, chain file (source URL, date accessed), input count, output count, unmapped count, unmapped file location.

#### Peak Filtering
Log: filter criteria (signalValue threshold, p-value cutoff), blacklist used and version, input/output region counts, what was removed.

#### Merging/Union Operations
Log: merge tool + version, merge distance parameter, input files (all accessions), sample tagging method, output count, overlap statistics.

#### R/Bioconductor Analysis
```json
{
  "tool": {
    "name": "DESeq2",
    "version": "1.42.0",
    "r_version": "4.3.2",
    "bioconductor_version": "3.18",
    "citation": "Love et al. 2014, Genome Biology, DOI:10.1186/s13059-014-0550-8"
  },
  "command": "DESeq2::results(dds, contrast=c('condition','treated','control'), alpha=0.05)",
  "parameters": {
    "design_formula": "~ batch + condition",
    "contrast": ["condition", "treated", "control"],
    "alpha": 0.05,
    "lfcThreshold": 0
  }
}
```

#### Python Analysis
```json
{
  "tool": {
    "name": "scanpy",
    "version": "1.9.6",
    "python_version": "3.11.5",
    "anndata_version": "0.10.3",
    "citation": "Wolf et al. 2018, Genome Biology, DOI:10.1186/s13059-017-1382-0"
  }
}
```

## Step 3: Record Software Environment

At the start of each analysis, capture the full environment:

### R Environment
```r
sessionInfo()
# Or more detailed:
devtools::session_info()
```

Log: R version, platform, attached packages with versions, loaded namespaces.

### Python Environment
```python
import pkg_resources
{pkg.key: pkg.version for pkg in pkg_resources.working_set}
```

Log: Python version, all installed packages with versions, virtual environment path.

### Command-Line Tools
For each tool used, record the version:
```bash
bedtools --version        # bedtools v2.31.0
samtools --version        # samtools 1.19
STAR --version            # 2.7.11a
macs2 --version           # macs2 2.2.9.1
liftOver                  # Kent tools (note: no --version flag; record binary date)
```

### System Information
```bash
uname -a                  # OS and kernel
nproc                     # CPU cores
free -h                   # Memory (Linux)
sysctl -n hw.memsize      # Memory (macOS)
nvidia-smi                # GPU info (if applicable)
```

## Step 4: Store Scripts

Every custom script used in the analysis should be stored in the `scripts/` directory with sequential numbering:

### Naming Convention
```
scripts/
├── 001_download_encode_data.sh
├── 002_filter_peaks.sh
├── 003_merge_samples.R
├── 004_chromhmm_segmentation.sh
├── 005_differential_analysis.R
└── 006_visualization.py
```

### Script Header Template
Every stored script should include a header:
```bash
#!/bin/bash
# Script: 002_filter_peaks.sh
# Project: H3K27ac analysis in human pancreas
# Date: 2024-01-15
# Author: Generated by ENCODE Connector
# Description: Filter H3K27ac peaks by signalValue and remove blacklisted regions
# Dependencies: bedtools v2.31.0, awk (GNU Awk 5.2.1)
# Input: ENCFF123ABC.bed (IDR thresholded peaks, GRCh38)
# Output: derived/filtered_peaks/H3K27ac_pancreas_filtered.bed
# Reference: hg38-blacklist.v2.bed.gz (Amemiya et al. 2019)
```

## Step 5: Log Derived Files to ENCODE Tracker

After each operation, register the derived file:
```
encode_log_derived_file(
    file_path="/path/to/output.bed",
    source_accessions=["ENCSR...", "ENCFF..."],
    description="H3K27ac peaks in pancreas, signalValue >= 4.5, blacklist-filtered",
    file_type="filtered_peaks",
    tool_used="bedtools v2.31.0 + awk",
    parameters="awk '$7 >= 4.5' | bedtools intersect -v blacklist"
)
```

Verify the provenance chain:
```
encode_get_provenance(file_path="/path/to/output.bed")
```

## Step 6: Version Control and Experiment Branching

### When the User Runs Multiple Versions
If the user tries different parameters or approaches:
1. Log EACH version as a separate operation with unique operation_id
2. Record what was different between versions
3. Ask the user which version to use going forward
4. Mark the chosen version as "selected" and others as "alternative"

### Example Version Log
```json
{
  "operation_id": "op_003a",
  "description": "Filter peaks - signalValue >= 4.5",
  "status": "alternative",
  "note": "Less stringent threshold, more peaks retained"
},
{
  "operation_id": "op_003b",
  "description": "Filter peaks - signalValue >= 7.0",
  "status": "selected",
  "note": "More stringent, user chose this for final analysis"
}
```

## Step 7: Auto-Generate Methods Sections

When the user requests methods writing, read the experiment log and generate publication-ready text.

### Methods Template Structure

**Data Acquisition**
> [Assay] data for [biosample] were obtained from the ENCODE Project (ENCODE Project Consortium 2020) via the ENCODE portal (https://www.encodeproject.org). [N] experiments were included (accessions: [list]). All experiments used [sequencer] with [read length] [SE/PE] reads, generating [N]M reads per replicate across [N] biological replicates. Data were processed by the ENCODE Uniform Processing Pipeline (version [X]).

**File Selection**
> [Output type] files aligned to [assembly] were selected for downstream analysis. Files were selected using ENCODE's preferred default designation. IDR thresholded peaks (Li et al. 2011) were used for [ChIP-seq/ATAC-seq] to ensure replicate concordance.

**Quality Assessment**
> Experiments were assessed for quality using ENCODE audit flags. Experiments with ERROR-level audits were excluded. ChIP-seq quality was evaluated using FRiP (≥[X]%), NSC (>[X]), RSC (>[X]), and NRF (≥[X]) metrics.

**Processing Steps**
For each operation in the log, generate a sentence:
> [Description]. [Tool] (version [X]; [citation]) was used with the following parameters: [parameters]. [Reference files] were obtained from [source] (version [X], accessed [date]). Of [N] input [regions/reads], [N] ([%]) passed filtering.

**Data Availability**
> All source data are available from the ENCODE portal under accessions [list]. Derived files, analysis scripts, and the complete provenance log are available at [repository URL]. Software versions: [list all tools and versions used].

### Scientific Documentation Standards

Methods sections MUST follow these principles for every computational step:

**Precision over approximation**
- Always use exact counts: "1,245 genes" not "~1,200 genes"; "42.3M reads" not "millions of reads"
- Every number that can be exact, should be exact

**Complete tool attribution**
- Every bioinformatics tool gets three things: name, version, and citation — no exceptions
- Don't just say "reads were aligned with STAR" — say which version, which index, which parameters

**Full reference specification**
- Always state BOTH genome build AND annotation version: "GRCh38.p14 with GENCODE v44" not just "hg38"
- Annotation versions matter: GENCODE v38 and v44 define different gene sets

**Experimental context**
- Name the sequencing platform: "Illumina HiSeq 4000" or "NovaSeq 6000"
- Report read characteristics: length, SE/PE, read count per sample, fragment size if relevant
- Report replicate details: number of biological replicates, sex breakdown (e.g., "3M/2F"), pooling strategy

**Show your filtering work**
- Every filtering step reports input count, pass count, fail count, and percentage: "Of 67,412 regions, 66,894 (99.2%) passed; 518 excluded"
- The reader should never wonder how much data was lost at any step

**Statistical rigor**
- Name the specific test AND the multiple testing correction: "Benjamini-Hochberg FDR < 0.05" not just "adjusted p-value"
- Name the software used for statistics with its version

**Data accessibility**
- Provide GEO/ENCODE accessions for both data deposition AND any external data reused
- Link to scripts, provenance logs, and derived files

**Orthogonal validation**
- Never rely on a single enrichment or pathway method — use 2+ complementary approaches so results don't depend on one algorithm's biases

### Citation Format for Tools
When generating methods, include proper citations:
| Tool | Citation |
|------|----------|
| bedtools | Quinlan & Hall 2010, Bioinformatics |
| samtools | Li et al. 2009, Bioinformatics |
| STAR | Dobin et al. 2013, Bioinformatics |
| featureCounts | Liao et al. 2014, Bioinformatics |
| edgeR | Robinson et al. 2010, Bioinformatics |
| MACS2 | Zhang et al. 2008, Genome Biology |
| DESeq2 | Love et al. 2014, Genome Biology |
| Seurat | Stuart et al. 2019, Cell |
| SCTransform | Hafemeister & Satija 2019, Genome Biology |
| CellRanger | 10x Genomics (cite version used) |
| Scanpy | Wolf et al. 2018, Genome Biology |
| ChromHMM | Ernst & Kellis 2012, Nature Methods |
| liftOver | Kent et al. 2002, Genome Research |
| HOMER | Heinz et al. 2010, Molecular Cell |
| deepTools | Ramirez et al. 2016, Nucleic Acids Research |
| Harmony | Korsunsky et al. 2019, Nature Methods |
| IDR | Li et al. 2011, Annals of Applied Statistics |
| WGCNA | Langfelder & Horvath 2008, BMC Bioinformatics |
| CibersortX | Newman et al. 2019, Nature Biotechnology |
| GSEA | Subramanian et al. 2005, PNAS |
| Gviz | Hahne & Ivanek 2016, Methods in Molecular Biology |
| GraphPad Prism | GraphPad Software (cite version) |
| DAVID | Huang et al. 2009, Nature Protocols |
| Enrichr | Kuleshov et al. 2016, Nucleic Acids Research |
| DEGAS | Li et al. 2022, Genome Biology |
| RRHO | Plaisier et al. 2010, Nucleic Acids Research |

## Step 8: Supplementary Data Tables

Generate supplementary tables following the scientific documentation standards above:

### Table S1: ENCODE Experiments Used
| Accession | Assay | Target | Biosample | Lab | Replicates | Sequencer | Read Length | Read Count | Library |
|-----------|-------|--------|-----------|-----|------------|-----------|-------------|------------|---------|
| ENCSR... | Histone ChIP-seq | H3K27ac | pancreas | Ren | 2 bio | HiSeq 4000 | 76bp PE | 42.3M | TruSeq |

### Table S2: Files Selected
| File Accession | Experiment | Format | Output Type | Assembly | Pipeline | Size | MD5 |
|---------------|-----------|--------|-------------|----------|----------|------|-----|
| ENCFF... | ENCSR... | bed narrowPeak | IDR thresholded peaks | GRCh38 | ENCODE v2.1 | 1.2MB | abc... |

### Table S3: Processing Steps
| Step | Description | Tool | Version | Input | Output | Parameters | Reference Files |
|------|------------|------|---------|-------|--------|------------|-----------------|
| 1 | Peak filtering | bedtools | 2.31.0 | ENCFF... | filtered.bed | signalValue≥4.5 | blacklist v2 |

### Table S4: Software Environment
| Software | Version | Citation |
|----------|---------|----------|
| R | 4.3.2 | R Core Team 2023 |
| Bioconductor | 3.18 | Huber et al. 2015 |
| DESeq2 | 1.42.0 | Love et al. 2014 |

Export using:
```
encode_export_data(format="csv")  # For Table S1
encode_get_citations(export_format="bibtex")  # For bibliography
```

## Pitfalls and Edge Cases

### Tool Version Drift
- Tool versions change over time; `bedtools v2.30` may produce different results than `v2.31`
- ALWAYS record versions at time of use, not at time of writing
- If re-running an analysis months later, verify tool versions match the log

### Reference File Versioning
- Genome annotations (GENCODE) release new versions regularly
- Chromosome size files differ between assemblies and even between UCSC/Ensembl conventions
- Blacklists have versions (v1 vs v2) that exclude different regions
- ALWAYS record the exact version and download URL/date

### Incomplete Provenance
- If a step was performed interactively (e.g., manual filtering in IGV), log it anyway with a note
- "No provenance" is worse than "approximate provenance"
- If a tool doesn't report its version, record the binary date and download source

### Multi-User Environments
- If multiple people work on the same project, log WHO performed each operation
- Use consistent file paths or relative paths in the log
- Store scripts in version control (git) alongside the provenance log

### Containerization for Exact Reproduction
- For maximum reproducibility, consider Docker/Singularity containers
- Record the container image and tag alongside the tool version
- Nextflow/Snakemake workflows can encode the full environment

## Walkthrough: Building a Complete Provenance Trail for an ENCODE Analysis

**Goal**: Document every step of an ENCODE analysis pipeline with full provenance — from raw data acquisition through processing, analysis, and derived outputs — enabling reproducibility and publication-ready methods.
**Context**: Reproducibility requires knowing exactly what data, tools, parameters, and versions produced each result. This skill automates provenance tracking.

### Step 1: Log data acquisition

```
encode_track_experiment(accession="ENCSR000AKA", notes="H3K27ac ChIP-seq in GM12878 for enhancer analysis")
```

Expected output:
```json
{
  "tracking": {
    "accession": "ENCSR000AKA",
    "action": "tracked"
  },
  "auto_linked_references": [
    {"type": "geo_accession", "id": "GSE29611"}
  ],
  "publications_found": 1,
  "publications": [
    {
      "pmid": "22955616",
      "doi": "10.1038/nature11247",
      "title": "An integrated encyclopedia of DNA elements in the human genome",
      "authors": "ENCODE Project Consortium",
      "journal": "Nature",
      "year": "2012",
      "abstract": ""
    }
  ],
  "pipelines_found": 0,
  "pipelines": []
}
```

The `notes` you pass are written to the tracker but not echoed here — read them back with `encode_list_tracked`.

### Step 2: Log file downloads with MD5 verification

```
encode_download_files(file_accessions=["ENCFF001ABC"], download_dir="/data/chipseq")
```

Expected output:
```json
{
  "downloaded": [
    {
      "accession": "ENCFF001ABC",
      "file_path": "/data/chipseq/ENCFF001ABC.bed.gz",
      "file_size": 1258291,
      "file_size_human": "1.2 MB",
      "success": true,
      "error": "",
      "md5_verified": true
    }
  ],
  "errors": [],
  "summary": {
    "total_requested": 1,
    "successful": 1,
    "failed": 0,
    "total_size": 1258291,
    "total_size_human": "1.2 MB"
  }
}
```

### Step 3: Log derived analysis outputs

```
encode_log_derived_file(
  source_accessions=["ENCFF001ABC", "ENCFF002DEF"],
  file_path="/data/analysis/gm12878_enhancers_filtered.bed",
  description="Filtered H3K27ac peaks: removed blacklist regions, merged within 500bp, filtered signalValue > 5",
  tool_used="bedtools v2.31.0",
  parameters="intersect -v (blacklist), merge -d 500, filter signalValue > 5"
)
```

Expected output:
```json
{
  "success": true,
  "record_id": 7,
  "file_path": "/data/analysis/gm12878_enhancers_filtered.bed",
  "source_accessions": ["ENCFF001ABC", "ENCFF002DEF"],
  "message": "Provenance logged. Use encode_get_provenance to view the full chain."
}
```

### Step 4: View the complete provenance chain

```
encode_get_provenance(file_path="/data/analysis/gm12878_enhancers_filtered.bed")
```

Expected output (the record also carries `created_at`, a float epoch timestamp):
```json
{
  "id": 7,
  "file_path": "/data/analysis/gm12878_enhancers_filtered.bed",
  "source_accessions": ["ENCFF001ABC", "ENCFF002DEF"],
  "description": "Filtered H3K27ac peaks: removed blacklist regions, merged within 500bp, filtered signalValue > 5",
  "file_type": "",
  "tool_used": "bedtools v2.31.0",
  "parameters": "intersect -v (blacklist), merge -d 500, filter signalValue > 5",
  "notes": "",
  "source_experiments": [
    {"accession": "ENCFF001ABC", "tracked": false},
    {"accession": "ENCFF002DEF", "tracked": false}
  ]
}
```

### Step 5: Generate provenance summary for publication

```
encode_summarize_collection()
```

**Interpretation**: The complete provenance chain enables automatic generation of methods sections: "H3K27ac ChIP-seq peaks (ENCFF001ABC) were filtered using ENCODE blacklist v2 (Amemiya et al. 2019) with bedtools v2.31.0..."

### Integration with downstream skills
- Provenance records from all skills feed into this skill for centralized tracking
- Provenance data supports → **scientific-writing** methods section generation
- File lineage connects to → **cite-encode** for proper ENCODE data attribution
- Pipeline provenance from → **pipeline-chipseq** through **pipeline-cutandrun** records processing steps

## Code Examples

### 1. Track an experiment for provenance
```
encode_track_experiment(accession="ENCSR000AKA", notes="GM12878 H3K27ac for enhancer catalog")
```

Expected output:
```json
{
  "tracking": {
    "accession": "ENCSR000AKA",
    "action": "tracked"
  },
  "publications_found": 0,
  "publications": [],
  "pipelines_found": 0,
  "pipelines": []
}
```

### 2. Log a derived analysis file
```
encode_log_derived_file(
  source_accessions=["ENCFF001ABC"],
  file_path="/data/peaks_filtered.bed",
  description="Blacklist-filtered peaks",
  tool_used="bedtools v2.31.0"
)
```

Expected output:
```json
{
  "success": true,
  "record_id": 8,
  "file_path": "/data/peaks_filtered.bed",
  "source_accessions": ["ENCFF001ABC"],
  "message": "Provenance logged. Use encode_get_provenance to view the full chain."
}
```

### 3. View provenance chain
```
encode_get_provenance(file_path="/data/peaks_filtered.bed")
```

Expected output:
```json
{
  "id": 8,
  "file_path": "/data/peaks_filtered.bed",
  "source_accessions": ["ENCFF001ABC"],
  "description": "Blacklist-filtered peaks",
  "file_type": "",
  "tool_used": "bedtools v2.31.0",
  "parameters": "",
  "notes": "",
  "source_experiments": [{"accession": "ENCFF001ABC", "tracked": false}]
}
```

## Integration

| This skill produces... | Feed into... | Using tool/skill |
|---|---|---|
| Provenance chain (accession → derived files) | Methods section generation | scientific-writing skill |
| Logged analysis steps with parameters | Reproducibility audit | publication-trust skill |
| MD5-verified file records | Data availability statement | cite-encode skill |
| Sequential script numbering | Pipeline documentation | pipeline-guide skill |
| Complete tool + version records | Tool citation list | cite-encode → BibTeX export |

## Related Skills

- `pipeline-guide` — ENCODE pipeline execution and monitoring
- `cite-encode` — Generating citations and bibliography for ENCODE data
- `quality-assessment` — Evaluating quality of ENCODE experiments
- `multi-omics-integration` — Multi-omics workflows that generate provenance
- `histone-aggregation` — Aggregation workflows that produce derived files
- `accessibility-aggregation` — ATAC-seq aggregation with provenance
- `geo-connector` — Log cross-references between ENCODE and GEO datasets
- `cross-reference` — Link experiments to PubMed, DOI, GEO, NCT IDs
- `publication-trust` — Verify literature claims backing analytical decisions

## Presenting Results

- Present provenance chain as: derived_file -> source_files -> ENCODE accessions, showing tools and parameters used. Include timestamps. Suggest: "Would you like to export this chain for your methods section?"

## For the request: "$ARGUMENTS"
