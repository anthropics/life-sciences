# Life Sciences Marketplace for Claude Code

This marketplace provides MCP (Model Context Protocol) servers and skills for life sciences tools. Install these plugins to access specialized research and analysis tools directly within Claude Code.

**What's included:**
- **MCP Servers**: Connect to external services like PubMed, BioRender, Synapse, and more
- **Skills**: Domain-specific workflows and analysis capabilities that extend Claude's expertise

## Quick Start

```bash
# Add the marketplace
/plugin marketplace add anthropics/life-sciences

# Install MCP servers
/plugin install pubmed@life-sciences
/plugin install biorender@life-sciences
/plugin install synapse@life-sciences
/plugin install wiley-scholar-gateway@life-sciences
/plugin install 10x-genomics@life-sciences
/plugin install encode-toolkit@life-sciences

# Install skills
/plugin install single-cell-rna-qc@life-sciences
/plugin install instrument-data-to-allotrope@life-sciences
/plugin install nextflow-development@life-sciences
/plugin install scvi-tools@life-sciences
/plugin install scientific-problem-selection@life-sciences
```

For servers requiring authentication (all except PubMed and ENCODE Toolkit), configure credentials after installation:
1. Type `/plugin` in Claude Code
2. Select "Manage plugins"
3. Find your installed server
4. Select "Configure"
5. Enter required credentials
6. Restart Claude Code

## Available Plugins

### Remote MCP Servers

#### PubMed
**Plugin ID**: `pubmed@life-sciences`

Search and access biomedical literature and research articles from PubMed.

**Requirements**: None - accessible to all users

#### BioRender
**Plugin ID**: `biorender@life-sciences`

Create and access scientific illustrations and diagrams.

**Requirements**: Free BioRender account (https://www.biorender.com)

#### Synapse.org
**Plugin ID**: `synapse@life-sciences`

Collaborative research data management platform by Sage Bionetworks.

**Requirements**: Free Synapse account (https://www.synapse.org)

#### Scholar Gateway (Wiley)
**Plugin ID**: `wiley-scholar-gateway@life-sciences`

Access academic research and publications from Wiley's Scholar Gateway.

**Requirements**: Free Scholar Gateway account

#### Consensus
**Plugin ID**: `consensus@life-sciences`

AI-powered search across 200M+ peer-reviewed scientific research papers, with evidence synthesis.

**Requirements**: Consensus account (https://consensus.app)

#### Cortellis Regulatory Intelligence (Clarivate)
**Plugin ID**: `cortellis@life-sciences`

Global drug regulatory intelligence covering submissions, approvals, and guidance documents.

**Requirements**: Cortellis subscription (https://clarivate.com/cortellis)

#### AdisInsight (Springer Nature)
**Plugin ID**: `adisinsight@life-sciences`

Drug development pipeline, clinical trials, safety, and deals intelligence.

**Requirements**: AdisInsight subscription (https://adisinsight.springer.com)

### Local MCP Servers (MCPB)

#### 10x Genomics Cloud
**Plugin ID**: `10x-genomics@life-sciences`

Access 10x Genomics Cloud analysis data and workflows.

**Requirements**:
- 10x Genomics Cloud account (https://www.10xgenomics.com/products/cloud-analysis)
- Access token (generate from: https://cloud.10xgenomics.com/account/security)
- Note: Only useful if you have analysis data in your account

### Local MCP Servers (stdio)

#### ENCODE Toolkit
**Plugin ID**: `encode-toolkit@life-sciences`

Search, download, track, and analyze functional genomics data from the ENCODE Project (https://www.encodeproject.org), the largest public catalog of functional genomic elements. Installs a local MCP server (launched with `uvx encode-toolkit==0.3.4`) together with 47 bundled skills.

**MCP server (20 tools):**
- Search experiments and files by assay (ChIP-seq, ATAC-seq, RNA-seq, WGBS, Hi-C, CUT&RUN, and more), organism, organ, biosample, and target
- Download files with automatic MD5 verification, with a dry-run preview for batch downloads
- Track experiments, citations (BibTeX/RIS), and derived-file provenance in a local SQLite database
- Link experiments to PubMed, bioRxiv, ClinicalTrials.gov, and GEO records for use alongside those MCP servers

**Bundled skills (47):**
- **Core**: setup, search, download, experiment tracking, cross-referencing
- **Analysis**: quality assessment, integrative analysis, regulatory elements, epigenome profiling, biosample comparison, visualization, motif analysis, peak annotation, batch analysis, functional screens (CRISPR, MPRA, STARR-seq)
- **Data aggregation**: histone marks, chromatin accessibility, Hi-C, DNA methylation
- **External databases**: GTEx, ClinVar, gnomAD, Ensembl, UCSC Genome Browser, GEO, GWAS Catalog, JASPAR, CELLxGENE
- **Workflows**: provenance, citations, variant annotation, single-cell, disease research, publication trust, tool installation, scientific writing, coordinate liftover, pipeline selection
- **Pipelines**: Nextflow + Docker pipelines for ChIP-seq, ATAC-seq, RNA-seq, WGBS, Hi-C, DNase-seq, and CUT&RUN
- **Meta-analysis**: scRNA-seq meta-analysis, multi-omics integration

**Requirements**:
- [uv](https://docs.astral.sh/uv/getting-started/installation/) installed locally (uv provisions Python 3.10+ automatically)
- Docker and Nextflow only if you run the pipeline skills
- No account needed for public ENCODE data; ENCODE access keys are optional and only used for unreleased or restricted datasets

**License**: The skills and manifest in this repository are Apache-2.0 (see `encode-toolkit/LICENSE.txt`). The MCP server is installed from PyPI at a pinned version and is AGPL-3.0; source at https://github.com/ammawla/encode-toolkit

### Skills

#### Single-Cell RNA-seq Quality Control
**Plugin ID**: `single-cell-rna-qc@life-sciences`

Automated quality control workflow for single-cell RNA-seq data following scverse best practices. Performs MAD-based filtering with comprehensive visualizations.

#### Instrument Data to Allotrope
**Plugin ID**: `instrument-data-to-allotrope@life-sciences`

Convert instrument data to Allotrope Simple Model (ASM) format for standardized data exchange and analysis.

#### Nextflow Development
**Plugin ID**: `nextflow-development@life-sciences`

Run nf-core bioinformatics pipelines (rnaseq, sarek, atacseq) on local or public GEO/SRA sequencing data. Designed for bench scientists who need to run large-scale omics analyses without specialized bioinformatics training.

**Supported pipelines:**
- **rnaseq**: Gene expression and differential expression analysis
- **sarek**: Germline and somatic variant calling (WGS/WES)
- **atacseq**: Chromatin accessibility analysis

**Features:**
- Download public datasets from GEO/SRA
- Auto-detect data types and suggest appropriate pipelines
- Generate pipeline-compatible samplesheets
- Environment validation and troubleshooting guidance

**Requirements**: Docker and Nextflow installed locally

#### scvi-tools
**Plugin ID**: `scvi-tools@life-sciences`

Deep learning toolkit for single-cell omics analysis using scvi-tools. Includes model selection guidance, training workflows, and integration pipelines for scVI, scANVI, totalVI, PeakVI, MultiVI, and more.

#### Scientific Problem Selection
**Plugin ID**: `scientific-problem-selection@life-sciences`

Systematic framework for scientific problem selection and strategic research decisions. Based on Fischbach & Walsh's methodology from Cell (2024), this skill helps researchers with project ideation, risk assessment, troubleshooting stuck projects, and strategic scientific planning.

**Use cases:**
- Pitch and refine new research ideas
- Evaluate project risks and feasibility
- Navigate decision trees in active projects
- Strategic research planning and problem choice

## Detailed Installation

### 1. Add the marketplace (one time)

```bash
/plugin marketplace add https://github.com/anthropics/life-sciences.git
```

### 2. Install specific plugins

```bash
# Remote MCP servers (no configuration needed for PubMed)
/plugin install pubmed@life-sciences
/plugin install biorender@life-sciences
/plugin install synapse@life-sciences
/plugin install wiley-scholar-gateway@life-sciences

# Local MCP servers (require configuration)
/plugin install 10x-genomics@life-sciences

# Local MCP servers (no configuration needed; requires uv)
/plugin install encode-toolkit@life-sciences

# Skills (no configuration needed)
/plugin install single-cell-rna-qc@life-sciences
/plugin install instrument-data-to-allotrope@life-sciences
/plugin install nextflow-development@life-sciences
/plugin install scvi-tools@life-sciences
/plugin install scientific-problem-selection@life-sciences
```

### 3. Configure credentials (if needed)

For servers requiring authentication, use the `/plugin` menu:
1. Type `/plugin` in Claude Code
2. Select "Manage plugins"
3. Find your installed server
4. Select "Configure" (if available)
5. Enter your API credentials

Or authenticate through the server's web interface when prompted.

### 4. Restart Claude Code

Restart to activate the MCP servers.

## Authentication Requirements

- **No authentication**: PubMed, ENCODE Toolkit
- **Free account required**: BioRender, Synapse, Wiley Scholar Gateway, Consensus
- **Paid/institutional account**: 10x Genomics (requires data in account to be useful), Cortellis, AdisInsight

## Support

For issues with:
- **Claude Code plugin system**: Report in #claude-cli-feedback on Anthropic Slack
- **Individual MCP servers**: Contact the respective provider's support

## License

Individual MCP servers are licensed by their respective providers. See each provider's terms of service for details.

## Removed Plugins

- **Benchling**: Removed because Benchling uses tenant-specific URLs which are not supported by the plugin system.
