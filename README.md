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
/plugin install biorxiv@life-sciences
/plugin install chembl@life-sciences
/plugin install clinical-trials@life-sciences
/plugin install open-targets@life-sciences
/plugin install owkin@life-sciences
/plugin install medidata@life-sciences
/plugin install tooluniverse@life-sciences

# Install skills
/plugin install single-cell-rna-qc@life-sciences
/plugin install instrument-data-to-allotrope@life-sciences
/plugin install nextflow-development@life-sciences
/plugin install scvi-tools@life-sciences
/plugin install clinical-trial-protocol@life-sciences
```

For servers requiring authentication (all except PubMed), configure credentials after installation:
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

#### bioRxiv/medRxiv
**Plugin ID**: `biorxiv@life-sciences`

Access to bioRxiv and medRxiv preprint data. Search and retrieve research papers in biological and medical sciences posted before peer review.

**Requirements**: None - accessible to all users

#### ChEMBL Database
**Plugin ID**: `chembl@life-sciences`

Access to the ChEMBL Database, a manually curated resource of bioactive drug-like compounds with quantitative binding and functional data against biological targets.

**Requirements**: None - accessible to all users

#### ClinicalTrials.gov
**Plugin ID**: `clinical-trials@life-sciences`

Access ClinicalTrials.gov data from the NIH/NLM registry of FDA-regulated clinical studies conducted worldwide.

**Requirements**: None - accessible to all users

#### Open Targets Platform
**Plugin ID**: `open-targets@life-sciences`

Drug target discovery and prioritisation platform that supports systematic identification of potential therapeutic drug targets. Integrates publicly available datasets to build and score target-disease associations.

**Requirements**: None - accessible to all users

#### Owkin
**Plugin ID**: `owkin@life-sciences`

Interact with AI agents built for biology to accelerate drug discovery and de-risk clinical trials. Currently powers HistoPLUS, an agent that transforms H&E slides from the TCGA database into granular, queryable insights for cell type quantification and spatial tumor microenvironment analysis.

**Requirements**: None - accessible to all users

#### Medidata
**Plugin ID**: `medidata@life-sciences`

Medidata's MCP Connector with Platform Help and Predictive Site Ranking for clinical trials. Platform Help provides access to Medidata platform documentation (Rave EDC, Data Connect, Clinical Data Studio). Predictive Site Ranking predicts which clinical trial sites align with enrollment goals during protocol planning.

**Requirements**: Medidata account for full functionality

### Local MCP Servers (MCPB)

#### 10x Genomics Cloud
**Plugin ID**: `10x-genomics@life-sciences`

Access 10x Genomics Cloud analysis data and workflows.

**Requirements**:
- 10x Genomics Cloud account (https://www.10xgenomics.com/products/cloud-analysis)
- Access token (generate from: https://cloud.10xgenomics.com/account/security)
- Note: Only useful if you have analysis data in your account

#### ToolUniverse
**Plugin ID**: `tooluniverse@life-sciences`

Democratizing AI scientists with ToolUniverse - AI agents and tools for scientific discovery and automation. Provides access to specialized tools for scientific research and computational analysis.

**Requirements**: None - local MCPB installation

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

#### Clinical Trial Protocol
**Plugin ID**: `clinical-trial-protocol@life-sciences`

Generate FDA/NIH-compliant clinical trial protocols for medical devices or drugs. Provides automated research of similar trials, FDA guidance documents, and comprehensive protocol generation following regulatory standards.

**Features:**
- Clinical trials research and FDA regulatory pathway analysis
- Protocol foundation, intervention, and operations sections
- Statistical sample size calculations
- NIH/FDA-compliant protocol templates
- Waypoint-based workflow for resumability

**Requirements**: Python dependencies (scipy, numpy) for statistical calculations

## Detailed Installation

### 1. Add the marketplace (one time)

```bash
/plugin marketplace add https://github.com/anthropics/life-sciences.git
```

### 2. Install specific plugins

```bash
# Remote MCP servers (no configuration needed for most)
/plugin install pubmed@life-sciences
/plugin install biorender@life-sciences
/plugin install synapse@life-sciences
/plugin install wiley-scholar-gateway@life-sciences
/plugin install biorxiv@life-sciences
/plugin install chembl@life-sciences
/plugin install clinical-trials@life-sciences
/plugin install open-targets@life-sciences
/plugin install owkin@life-sciences
/plugin install medidata@life-sciences

# Local MCP servers (require configuration)
/plugin install 10x-genomics@life-sciences
/plugin install tooluniverse@life-sciences

# Skills (no configuration needed except clinical-trial-protocol)
/plugin install single-cell-rna-qc@life-sciences
/plugin install instrument-data-to-allotrope@life-sciences
/plugin install nextflow-development@life-sciences
/plugin install scvi-tools@life-sciences
/plugin install clinical-trial-protocol@life-sciences
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

- **No authentication**: PubMed, bioRxiv, ChEMBL, ClinicalTrials.gov, Open Targets, Owkin, ToolUniverse
- **Free account required**: BioRender, Synapse, Wiley Scholar Gateway
- **Paid/institutional account**: 10x Genomics (requires data in account to be useful), Medidata (requires Medidata account)

## Support

For issues with:
- **Claude Code plugin system**: Report in #claude-cli-feedback on Anthropic Slack
- **Individual MCP servers**: Contact the respective provider's support

## License

Individual MCP servers are licensed by their respective providers. See each provider's terms of service for details.

## Removed Plugins

- **Benchling**: Removed because Benchling uses tenant-specific URLs which are not supported by the plugin system.
