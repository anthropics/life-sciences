# Life Sciences Marketplace for Claude Code

This marketplace provides skills for life sciences research. Skills are domain-specific workflows that extend Claude's expertise.

## Quick Start

```bash
# Add the marketplace
/plugin marketplace add anthropics/life-sciences

# Install skills
/plugin install clinical-trial-protocol@life-sciences
/plugin install scientific-problem-selection@life-sciences
/plugin install binder-gauntlet@life-sciences
```

## Available Skills

### Clinical Trial Protocol
**Plugin ID**: `clinical-trial-protocol@life-sciences`

Generate FDA/NIH-compliant clinical trial protocols for medical devices or drugs.

### Scientific Problem Selection
**Plugin ID**: `scientific-problem-selection@life-sciences`

Systematic framework for scientific problem selection and strategic research decisions. Based on Fischbach & Walsh's methodology from Cell (2024), this skill helps researchers with project ideation, risk assessment, troubleshooting stuck projects, and strategic scientific planning.

**Use cases:**
- Pitch and refine new research ideas
- Evaluate project risks and feasibility
- Navigate decision trees in active projects
- Strategic research planning and problem choice

### Clinical Trial Forecast
**Plugin ID**: `clinical-trial-forecast@life-sciences`

Estimate the probability that a registered clinical trial meets its goal, with a cited report and the simulation model behind the number. Runs in Claude Science; see the [skill's README](clinical-trial-forecast/README.md).

### Binder Gauntlet
**Plugin ID**: `binder-gauntlet@life-sciences`

Validate de novo designed protein binders with independent re-prediction and controls, profile their developability, and write a dossier that states what was not checked. All results are computational predictions.

## Removed Plugins

This repository has been refocused on skills, so it no longer carries MCP connector plugins. If you installed one, it will no longer update from this marketplace. Earlier versions remain available in the git history.

The following skills were also removed: `instrument-data-to-allotrope`, `nextflow-development`, `scvi-tools` and `single-cell-rna-qc`.
