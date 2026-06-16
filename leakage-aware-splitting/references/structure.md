# Structural splitting (Foldseek)

Use when the task is structure/fold/interface-centric, or when remote homologs that
share a fold but not detectable sequence identity are expected. A sequence-identity
split leaks in these cases: two proteins below 30% identity can be structural
near-duplicates (PDB 1K3F and 1K9S share only 26.5% identity yet are near-identical
structurally), so an interaction or fold label leaks across a sequence split.

## When to use

- Docking, interface prediction, fold classification, binding-pocket tasks.
- You need to test generalization to **new folds** (CASP "free modeling").
- Remote homology with low sequence identity is plausible in your data.

## Recipe

```bash
foldseek easy-cluster structures/ clu tmp \
    -c 0.8 --tmscore-threshold 0.5
```

- **TM-score >= 0.5** is the standard fold-equivalence threshold; Foldseek also uses
  E-value < 0.01. Cluster, then assign whole structural clusters to one split.
- Foldseek encodes structures in a 3Di alphabet and searches via MMseqs2, reaching
  TM-align/DALI sensitivity thousands of times faster, so all-vs-all structural
  leakage screening is tractable.

`scripts/cluster_structure.py:foldseek_cluster` wraps this and returns per-item
cluster labels. There is **no pure-Python fallback** (real structural similarity
needs 3D coordinates) - if Foldseek is absent, either supply precomputed cluster
labels or fall back to a sequence split and disclose that structural leakage is
unchecked.

## Diagnostics

- For each test structure, the max TM-score (or min Foldseek E-value) to any train
  structure; assert below threshold.
- For protein-protein interfaces, split on **interface** structural similarity
  (iDist/iAlign), not just monomer structure - DIPS suffers ~53% structural leakage
  even under sequence-family splitting.

## When sequence is enough

If the label tracks sequence (signal peptides, many sequence-only property tasks) and
the data has no structures, a sequence split suffices. Add structure only when the
task or data demands it.

## Citations

- van Kempen et al., Foldseek, Nat. Biotechnol. 2023.
- Barrio-Hernandez et al., Foldseek cluster of the AlphaFold DB, Nature 622:637, 2023.
- Bushuiev et al., PPI leakage, arXiv 2404.10457, 2024 (iDist; DIPS 53% leakage).
