# Life Sciences Marketplace for Claude Code

This marketplace provides MCP (Model Context Protocol) servers and skills for life sciences tools. Install these plugins to access specialized research and analysis tools directly within Claude Code.

**What's included:**
- **MCP Servers**: Connect to external services like PubMed, BioRender, Benchling, and more
- **Skills**: Domain-specific workflows and analysis capabilities that extend Claude's expertise

## Quick Start

```bash
# Add the marketplace
claude marketplace add https://github.com/anthropics/life-sciences.git

# Install MCP servers
claude plugin install pubmed@life-sciences
claude plugin install biorender@life-sciences
claude plugin install synapse@life-sciences
claude plugin install wiley-scholar-gateway@life-sciences
claude plugin install benchling-mcp@life-sciences
claude plugin install 10x-genomics@life-sciences

# Install skills
claude plugin install single-cell-rna-qc@life-sciences
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

### Local MCP Servers (MCPB)

#### Benchling
**Plugin ID**: `benchling-mcp@life-sciences`

Access Benchling notebooks, entries, schemas, and more.

**Requirements**:
- Benchling account with Benchling AI, API access, and "Ask" functionality enabled
- API key (generate from Benchling Settings)
- Benchling tenant URL

#### 10x Genomics Cloud
**Plugin ID**: `10x-genomics@life-sciences`

Access 10x Genomics Cloud analysis data and workflows.

**Requirements**:
- 10x Genomics Cloud account (https://www.10xgenomics.com/products/cloud-analysis)
- Access token (generate from: https://cloud.10xgenomics.com/account/security)
- Note: Only useful if you have analysis data in your account

### Skills

#### Single-Cell RNA-seq Quality Control
**Plugin ID**: `single-cell-rna-qc@life-sciences`

Automated quality control workflow for single-cell RNA-seq data following scverse best practices. Performs MAD-based filtering with comprehensive visualizations.

**Requirements**:
- Python packages: anndata, scanpy, scipy, matplotlib, seaborn, numpy
- Input data: `.h5ad` (AnnData) or `.h5` (10X Genomics) files

**Installation**:
```bash
claude plugin install single-cell-rna-qc@life-sciences
```

## Detailed Installation

### 1. Add the marketplace (one time)

```bash
claude marketplace add https://github.com/anthropics/life-sciences.git
```

### 2. Install specific plugins

```bash
# Remote MCP servers (no configuration needed for PubMed)
claude plugin install pubmed@life-sciences
claude plugin install biorender@life-sciences
claude plugin install synapse@life-sciences
claude plugin install wiley-scholar-gateway@life-sciences

# Local MCP servers (require configuration)
claude plugin install benchling-mcp@life-sciences
claude plugin install 10x-genomics@life-sciences

# Skills (no configuration needed)
claude plugin install single-cell-rna-qc@life-sciences
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

- **No authentication**: PubMed
- **Free account required**: BioRender, Synapse, Wiley Scholar Gateway
- **Paid/institutional account**: Benchling (requires tenant with specific features), 10x Genomics (requires data in account to be useful)

## Support

For issues with:
- **Claude Code plugin system**: Report in #claude-cli-feedback on Anthropic Slack
- **Individual MCP servers**: Contact the respective provider's support

## License

Individual MCP servers are licensed by their respective providers. See each provider's terms of service for details.
