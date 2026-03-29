# CRISPR OS

Evidence-first CRISPR experiment design skill submitted from [`Hordago-Labs/crispr-os`](https://github.com/Hordago-Labs/crispr-os).

## Overview

This skill packages the marketplace-facing entry point, reference notes, and worked examples for CRISPR OS. The upstream project provides the full deterministic runtime for:

- guide discovery and ranking
- off-target review
- primer and oligo design
- validation planning
- provenance-tracked experiment reporting

## Included Here

- `SKILL.md` for skill routing and usage guidance
- `references/overview.md` for capability summary and quick start
- `references/evidence-first-architecture.md` for the deterministic design pattern
- `references/workflow-examples.md` for adapted knockout, base editing, and CRISPRa examples

## Upstream Runtime

Run the full system from the source repository:

```bash
git clone https://github.com/Hordago-Labs/crispr-os.git
cd crispr-os
python3 scripts/experiment_runner.py run "Knock out TP53 in HEK293T cells"
```
