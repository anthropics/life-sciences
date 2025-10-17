# Life Sciences MCP Servers for Claude Code

This marketplace provides MCP (Model Context Protocol) servers for life sciences tools. Install these plugins to access specialized research and analysis tools directly within Claude Code.

## Available Servers

### Remote MCP Servers

#### PubMed
**Plugin**: `pubmed@life-sciences`

Search and access biomedical literature and research articles from PubMed.

**Requirements**: None - accessible to all users

#### BioRender
**Plugin**: `biorender@life-sciences`

Create and access scientific illustrations and diagrams.

**Requirements**: Free BioRender account (https://www.biorender.com)

#### Synapse.org
**Plugin**: `synapse@life-sciences`

Collaborative research data management platform by Sage Bionetworks.

**Requirements**: Free Synapse account (https://www.synapse.org)

#### Scholar Gateway (Wiley)
**Plugin**: `wiley-scholar-gateway@life-sciences`

Access academic research and publications from Wiley's Scholar Gateway.

**Requirements**: Free Scholar Gateway account

### Local MCP Servers (MCPB)

#### Benchling
**Plugin**: `benchling-mcp@life-sciences`

Access Benchling notebooks, entries, schemas, and more.

**Requirements**:
- Benchling account with Benchling AI, API access, and "Ask" functionality enabled
- API key (generate from Benchling Settings)
- Benchling tenant URL

#### 10x Genomics Cloud
**Plugin**: `10x-genomics@life-sciences`

Access 10x Genomics Cloud analysis data and workflows.

**Requirements**:
- 10x Genomics Cloud account (https://www.10xgenomics.com/products/cloud-analysis)
- Access token (generate from: https://cloud.10xgenomics.com/account/security)
- Note: Only useful if you have analysis data in your account

## Installation

### Quick Start - Install All Servers

```bash
claude marketplace add https://github.com/anthropics/life-sciences.git
claude plugin install life-sciences pubmed
claude plugin install life-sciences biorender
claude plugin install life-sciences synapse
claude plugin install life-sciences wiley-scholar-gateway
claude plugin install life-sciences benchling-mcp
claude plugin install life-sciences 10x-genomics
```

### Individual Installation

1. **Add the marketplace** (one time):
```bash
claude marketplace add https://github.com/anthropics/life-sciences.git
```

2. **Install specific servers**:
```bash
# Remote servers (no configuration needed for PubMed)
claude plugin install life-sciences pubmed
claude plugin install life-sciences biorender
claude plugin install life-sciences synapse
claude plugin install life-sciences wiley-scholar-gateway

# Local servers (require configuration)
claude plugin install life-sciences benchling-mcp
claude plugin install life-sciences 10x-genomics
```

3. **Configure credentials** (for servers requiring authentication):

After installation, configure your credentials via the `/plugin` menu:
1. Type `/plugin` in Claude Code
2. Select "Manage plugins"
3. Find your installed server
4. Select "Configure" (if available)
5. Enter your API credentials

Or authenticate through the server's web interface when prompted.

4. **Restart Claude Code** to activate the MCP servers

## Server Details

### Authentication Requirements

- **No authentication**: PubMed
- **Free account required**: BioRender, Synapse, Wiley Scholar Gateway
- **Paid/institutional account**: Benchling (requires tenant with specific features), 10x Genomics (requires data in account to be useful)

## Support

For issues with:
- **Claude Code plugin system**: Report in #claude-cli-feedback on Anthropic Slack
- **Individual MCP servers**: Contact the respective provider's support

## License

Individual MCP servers are licensed by their respective providers. See each provider's terms of service for details.
