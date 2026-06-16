#!/usr/bin/env python3
"""Protein-ligand (drug-target) splitting: control leakage on BOTH axes.

Under a random split, deep DTI models reach AUROC > 0.98 - but this reflects hidden
ligand bias (correct predictions from drug features alone, not interaction patterns).
The honest split is "double-cold" (cold-drug + cold-target): test pairs whose protein
AND ligand are both unseen in training. PDBbind->CASF leakage is the canonical case
(>700 training complexes share near-duplicate similarity with CASF, ~45% of test
complexes), and de-leaking it drops top models toward baseline.

This module implements a greedy double-cold split as a dependency-free fallback and
documents DataSAIL as the production tool that solves the 2D assignment optimally.

A non-reported protein-ligand pair is *untested*, not a confirmed negative - keep
that in mind when constructing negatives; this splitter only partitions the pairs you
give it.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Sequence


def double_cold_split(
    protein_cluster: Sequence[int],
    ligand_cluster: Sequence[int],
    fractions: dict[str, float] | None = None,
    seed: int = 0,
) -> tuple[list[str | None], dict]:
    """Greedy cold-drug + cold-target split.

    Each pair i has a protein cluster and a ligand cluster. We partition the set of
    protein clusters and (independently) the set of ligand clusters into the requested
    splits. A pair is assigned to split S only if BOTH its protein cluster and its
    ligand cluster were assigned to S; pairs falling in the off-diagonal blocks
    (protein in train, ligand in test, etc.) are DISCARDED, because keeping them would
    leak one axis. Discarding is inherent to double-cold splits - the function reports
    how much was dropped so the caller can disclose it.

    Returns (assignment, info) where assignment[i] is a split name or None (discarded).
    """
    if fractions is None:
        fractions = {"train": 0.8, "test": 0.2}
    fractions = {k: v for k, v in fractions.items() if v > 0}
    splits = list(fractions.keys())

    prot_clusters = sorted(set(protein_cluster))
    lig_clusters = sorted(set(ligand_cluster))

    prot_assign = _assign_groups_by_fraction(prot_clusters, fractions, seed)
    lig_assign = _assign_groups_by_fraction(lig_clusters, fractions, seed + 1)

    assignment: list[str | None] = []
    kept = {s: 0 for s in splits}
    discarded = 0
    for p, l in zip(protein_cluster, ligand_cluster):
        ps, ls = prot_assign[p], lig_assign[l]
        if ps == ls:
            assignment.append(ps)
            kept[ps] += 1
        else:
            assignment.append(None)
            discarded += 1

    n = len(protein_cluster)
    info = {
        "kept_per_split": kept,
        "discarded": discarded,
        "discarded_fraction": discarded / n if n else 0.0,
        "n_protein_clusters": len(prot_clusters),
        "n_ligand_clusters": len(lig_clusters),
        "note": "off-diagonal pairs discarded to keep both axes cold; "
                "use DataSAIL for an optimal 2D assignment that retains more data",
    }
    return assignment, info


def _assign_groups_by_fraction(clusters: Sequence[int], fractions: dict[str, float], seed: int) -> dict[int, str]:
    """Deterministically assign cluster ids to splits hitting the fractions by count."""
    splits = list(fractions.keys())
    targets = {s: fractions[s] * len(clusters) for s in splits}
    counts = {s: 0 for s in splits}
    out: dict[int, str] = {}
    # Deterministic order, lightly perturbed by seed without a global RNG.
    ordered = sorted(clusters, key=lambda c: (hash((seed, c)) & 0xFFFFFFFF, c))
    for c in ordered:
        best = max(splits, key=lambda s: (targets[s] - counts[s], -splits.index(s)))
        out[c] = best
        counts[best] += 1
    return out


def leakage_across_axes(
    protein_cluster: Sequence[int],
    ligand_cluster: Sequence[int],
    assignment: Sequence[str | None],
) -> dict:
    """Verify double-cold integrity: no protein or ligand cluster shared train<->test."""
    prot_to_splits = defaultdict(set)
    lig_to_splits = defaultdict(set)
    for p, l, s in zip(protein_cluster, ligand_cluster, assignment):
        if s is None:
            continue
        prot_to_splits[p].add(s)
        lig_to_splits[l].add(s)
    prot_leaks = {p: sp for p, sp in prot_to_splits.items() if len(sp) > 1}
    lig_leaks = {l: sp for l, sp in lig_to_splits.items() if len(sp) > 1}
    return {
        "protein_axis_clean": not prot_leaks,
        "ligand_axis_clean": not lig_leaks,
        "protein_leaks": len(prot_leaks),
        "ligand_leaks": len(lig_leaks),
    }
