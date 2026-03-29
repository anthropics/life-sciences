# Evidence-First Architecture

CRISPR OS is organized around a constrained-LLM pattern:

1. Claude parses researcher intent into a structured experiment description.
2. Deterministic scripts and scientific engines perform the computation.
3. Validation gates verify the resulting artifacts.
4. Provenance manifests record hashes, parameters, and tool versions.

## Why This Matters

- Errors should fail visibly, not appear as plausible prose.
- Identical inputs should produce identical outputs whenever possible.
- Scientific recommendations should point back to concrete artifacts and scoring rules.

## Representative Upstream Claims

Adapted from the upstream project documentation:

- The system treats outputs as compiled scientific artifacts rather than free-form AI text.
- Guide design, plate layouts, validation planning, and provenance are part of a deterministic runtime.
- Claude is used at decision boundaries, not as the primary scientific calculator.

## Practical Consequence

When using this skill in Claude Code, route toward explicit steps such as:

- select modality
- review safety constraints
- run guide-design scripts
- check off-target results
- generate validation plan
- package manifest and report
