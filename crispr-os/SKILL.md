---
name: crispr-os
description: Evidence-first CRISPR experiment design. Use when users need (1) guide RNA design for knockout, base editing, CRISPRa/CRISPRi, or prime editing, (2) off-target review and safety checks, (3) primer or oligo design, (4) validation planning across DNA/RNA/protein/phenotype tiers, (5) plate layout and protocol generation, or (6) provenance-tracked experiment planning and reporting. Triggers include CRISPR, sgRNA, gRNA, guide RNA, knockout, base editing, prime editing, CRISPRa, CRISPRi, PAM, SpCas9, Cas12a, off-target, validation assay, genotyping primer, and experiment manifest.
---

# CRISPR OS

CRISPR OS is an evidence-first CRISPR experiment design skill adapted from [`Hordago-Labs/crispr-os`](https://github.com/Hordago-Labs/crispr-os). It is designed around a deterministic scientific runtime: Claude handles intent and recovery, while guide finding, GC/Tm checks, off-target scoring, validation planning, and provenance generation are performed by scripts and external scientific engines.

## What This Skill Covers

- Knockout guide design with PAM scanning, GC/Tm checks, and ranking
- Base editing, CRISPRa/CRISPRi, and prime editing workflow selection
- Off-target review, delivery planning, and safety/biosafety checks
- Validation planning across DNA, RNA, protein, and phenotype tiers
- Experiment packaging with manifests, reports, and worked examples

## Important Scope Note

This marketplace package ships the skill entry, references, and usage examples. The full executable runtime lives in the upstream repository:

- Repository: `https://github.com/Hordago-Labs/crispr-os`
- Authoritative runtime: `scripts/experiment_runner.py`
- Primary skill source: root `SKILL.md` in the upstream repo

## Quick Start

```bash
git clone https://github.com/Hordago-Labs/crispr-os.git
cd crispr-os

# Design-only mode
python3 scripts/experiment_runner.py run "Knock out TP53 in HEK293T cells"

# Full mode with analysis and report generation
python3 scripts/experiment_runner.py run "Knock out TP53 in HEK293T cells" --mode full
```

For the exact workflow patterns and examples packaged with this marketplace submission, read:

- `references/overview.md`
- `references/evidence-first-architecture.md`
- `references/workflow-examples.md`

## When to Use This Skill

Use this skill when users ask for:

- CRISPR guide design or ranking
- Knockout, base editing, CRISPRa, CRISPRi, or prime editing design help
- Off-target assessment and guide safety review
- Validation assay planning for edited cells
- Provenance-tracked CRISPR experiment planning

## Core Expectations

1. Read the workflow references before choosing a modality.
2. Apply biosafety review before producing guide suggestions.
3. Prefer deterministic outputs and explicit scoring criteria over free-form recommendations.
4. Treat the upstream runtime as authoritative for any executable design or report workflow.
