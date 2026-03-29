# CRISPR OS Overview

CRISPR OS is an evidence-first CRISPR design system built for deterministic experiment planning. The upstream project combines stdlib-only helper scripts with external scientific engines so that scoring, filtering, and manifest generation are reproducible and auditable.

## Primary Capabilities

- Guide RNA design for knockout workflows
- Modality routing for base editing, CRISPRa/CRISPRi, and prime editing
- GC content, homopolymer, and melting temperature checks
- Off-target review and safety-oriented filtering
- Validation planning across DNA, RNA, protein, and phenotype tiers
- Provenance manifests and report packaging

## Authoritative Runtime

The upstream project exposes the end-to-end workflow via:

```bash
python3 scripts/experiment_runner.py run "Knock out TP53 in HEK293T cells"
python3 scripts/experiment_runner.py run "Knock out TP53 in HEK293T cells" --mode full
```

The design pipeline referenced by the source skill includes:

- `scripts/pam_scanner.py`
- `scripts/gc_content.py`
- `scripts/tm_calculator.py`
- `scripts/off_target_scorer.py`
- `scripts/validation_planner.py`
- `scripts/provenance.py`

## Safety Expectations

- Perform biosafety review before proposing guides
- Block germline editing or dual-use requests
- Include scoring methods and search parameters with every recommendation
- Prefer loud validation failures over speculative output
