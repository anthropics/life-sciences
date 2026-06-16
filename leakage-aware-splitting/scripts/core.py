#!/usr/bin/env python3
"""Modality-agnostic, pure-Python core for leakage-aware splitting.

The scientifically important part of a leakage-aware split is NOT the clustering
tool (MMseqs2, Foldseek, RDKit, ...) but the rule that whole *clusters* of
similar items are assigned to exactly one split. This module implements that
rule independently of any external tool, so it is deterministic, hermetic, and
testable without binaries.

Pipeline:
    similar pairs (from any modality)  ->  clusters (single-linkage / union-find)
    clusters + target fractions        ->  split assignment (whole clusters kept together)

Tool-specific wrappers (cluster_sequence.py, cluster_structure.py, ...) only need
to produce either (a) a list of similar pairs above a threshold, or (b) explicit
cluster labels. Everything downstream is shared and lives here.

Why single-linkage (connected components) for the default clustering: with
centroid-based clustering (CD-HIT / MMseqs2 set-cover) two sequences placed in
different clusters can still be similar to each other (each was closer to its own
representative), so assigning whole clusters to splits does NOT by itself
guarantee between-split dissimilarity. Connected-component clustering uses
transitive reachability (A~B, B~C => {A,B,C}), giving the strongest guarantee
that members of different clusters fall below threshold.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Sequence

# Canonical split names. Three-way is the default; two-way drops "val".
SPLIT_NAMES = ("train", "val", "test")


class UnionFind:
    """Disjoint-set with path compression and union by rank."""

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        # path compression
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def clusters_from_pairs(n_items: int, similar_pairs: Iterable[tuple[int, int]]) -> list[int]:
    """Single-linkage (connected-component) clustering from similar pairs.

    Parameters
    ----------
    n_items : int
        Total number of items (indices 0..n_items-1). Items with no edges are singletons.
    similar_pairs : iterable of (i, j)
        Index pairs whose similarity is at/above the leakage threshold. Order and
        direction do not matter; duplicates are harmless.

    Returns
    -------
    list[int]
        cluster_labels[i] is the cluster id of item i. Labels are contiguous
        integers 0..n_clusters-1, assigned in order of first appearance for
        determinism.
    """
    uf = UnionFind(n_items)
    for i, j in similar_pairs:
        uf.union(i, j)
    # Renumber roots to contiguous, deterministic labels.
    label_of_root: dict[int, int] = {}
    labels = [0] * n_items
    for i in range(n_items):
        root = uf.find(i)
        if root not in label_of_root:
            label_of_root[root] = len(label_of_root)
        labels[i] = label_of_root[root]
    return labels


@dataclass
class SplitResult:
    """A complete split assignment plus the cluster structure it respects."""

    assignment: list[str]  # assignment[i] in SPLIT_NAMES for item i
    cluster_labels: list[int]
    fractions: dict[str, float]
    realized_fractions: dict[str, float] = field(default_factory=dict)
    n_clusters: int = 0
    seed: int = 0

    def indices(self, split: str) -> list[int]:
        return [i for i, s in enumerate(self.assignment) if s == split]


def assign_clusters_to_splits(
    cluster_labels: Sequence[int],
    fractions: dict[str, float] | None = None,
    seed: int = 0,
    stratify_labels: Sequence | None = None,
) -> SplitResult:
    """Assign whole clusters to splits, hitting target fractions as closely as possible.

    Guarantees by construction: every item lands in exactly one split, and no
    cluster id appears in more than one split (cluster integrity). Because cluster
    sizes vary, the realized item-level fractions will differ from the requested
    ones; this is expected and reported, not an error.

    Algorithm: a deterministic greedy bin-packing in two steps. First, when there
    are at least as many clusters as splits, the smaller splits (all but the
    largest-target one) are *seeded* with the smallest available clusters so they
    cannot be starved into emptiness; with whole-cluster assignment an empty split is
    otherwise a real failure mode (a few large clusters and the small splits never get
    one). Seeding with the smallest clusters guarantees each small split is non-empty
    while barely perturbing the ratio. Second, the remaining clusters are sorted
    largest-first (ties by cluster id) and each is placed in the split currently most
    under its target count; placing large clusters first keeps them from overshooting a
    small split late. The seed only perturbs tie handling among equal-size clusters via
    a stable shuffle, keeping runs with the same seed byte-identical. When there are
    fewer clusters than splits the request is unsatisfiable and some split is
    unavoidably empty; that is left for the no_empty_splits invariant to report rather
    than silently rebalanced (which would break cluster integrity).

    Parameters
    ----------
    cluster_labels : sequence of int
        Per-item cluster id (e.g. from clusters_from_pairs).
    fractions : dict, optional
        Target item-level fractions per split. Defaults to 80/10/10. Must sum to 1
        (within tolerance). Splits with fraction 0 are dropped.
    seed : int
        Reproducibility seed.
    stratify_labels : sequence, optional
        Per-item class labels. When given, clusters are assigned within strata
        defined by each cluster's majority label, so the class balance is
        preserved across splits as far as cluster granularity allows.

    Returns
    -------
    SplitResult
    """
    n = len(cluster_labels)
    if fractions is None:
        fractions = {"train": 0.8, "val": 0.1, "test": 0.1}
    fractions = {k: v for k, v in fractions.items() if v > 0}
    total = sum(fractions.values())
    if abs(total - 1.0) > 1e-6:
        raise ValueError(f"fractions must sum to 1.0, got {total}")

    # Group item indices by cluster.
    members: dict[int, list[int]] = defaultdict(list)
    for idx, c in enumerate(cluster_labels):
        members[c].append(idx)

    # Determine a stratum key per cluster (majority class label, or a single shared
    # stratum when not stratifying).
    if stratify_labels is not None:
        if len(stratify_labels) != n:
            raise ValueError("stratify_labels must align with cluster_labels")
        cluster_stratum: dict[int, object] = {}
        for c, idxs in members.items():
            counts: dict[object, int] = defaultdict(int)
            for i in idxs:
                counts[stratify_labels[i]] += 1
            # Majority label; ties broken by sorted label for determinism.
            cluster_stratum[c] = max(sorted(counts), key=lambda k: counts[k])
    else:
        cluster_stratum = {c: None for c in members}

    split_names = list(fractions.keys())
    assignment = [None] * n  # type: ignore[var-annotated]
    counts = {s: 0 for s in split_names}

    # Process strata in deterministic order so stratified runs are reproducible.
    strata = sorted({cluster_stratum[c] for c in members}, key=lambda x: (x is None, str(x)))
    for stratum in strata:
        stratum_clusters = [c for c in members if cluster_stratum[c] == stratum]
        n_stratum = sum(len(members[c]) for c in stratum_clusters)
        # Per-stratum targets so each class is split in the requested proportions.
        targets = {s: fractions[s] * n_stratum for s in split_names}
        local_counts = {s: 0 for s in split_names}

        # Step 1: seed the smaller splits (all but the largest-target one) with the
        # smallest clusters so they are never starved into emptiness when avoidable.
        # Requires at least one cluster per split; otherwise some split is unavoidably
        # empty and the no_empty_splits invariant reports it.
        seeded: set = set()
        if len(split_names) > 1 and len(stratum_clusters) >= len(split_names):
            by_target = sorted(split_names, key=lambda s: (targets[s], split_names.index(s)))
            smallest_first = sorted(stratum_clusters, key=lambda c: (len(members[c]), c))
            for s, c in zip(by_target[:-1], smallest_first):  # exclude the largest-target split
                for i in members[c]:
                    assignment[i] = s
                local_counts[s] += len(members[c])
                counts[s] += len(members[c])
                seeded.add(c)

        # Step 2: largest-first deficit fill for the rest. Sort largest-first with a
        # deterministic tie-break by cluster id; seed enters only as a stable rotation
        # of equal-size groups, preserving determinism.
        ordered = sorted(stratum_clusters, key=lambda c: (-len(members[c]), c))
        if seed:
            ordered = _stable_seeded_reorder_within_ties(ordered, members, seed)

        for c in ordered:
            if c in seeded:
                continue
            size = len(members[c])
            # Choose the split with the largest remaining deficit (target - current),
            # tie-broken by the split's global order for reproducibility.
            best = max(
                split_names,
                key=lambda s: (targets[s] - local_counts[s], -split_names.index(s)),
            )
            for i in members[c]:
                assignment[i] = best
            local_counts[best] += size
            counts[best] += size

    realized = {s: counts[s] / n if n else 0.0 for s in split_names}
    return SplitResult(
        assignment=assignment,  # type: ignore[arg-type]
        cluster_labels=list(cluster_labels),
        fractions=fractions,
        realized_fractions=realized,
        n_clusters=len(members),
        seed=seed,
    )


def _stable_seeded_reorder_within_ties(ordered, members, seed):
    """Deterministically reorder only equal-size runs, keyed by seed.

    Keeps the largest-first ordering intact (which drives quality) while letting
    the seed vary which of two equal-size clusters is placed first. Uses a hash of
    (seed, cluster_id) so the result is reproducible for a given seed and never
    calls a global RNG.
    """
    out = []
    i = 0
    while i < len(ordered):
        j = i
        while j < len(ordered) and len(members[ordered[j]]) == len(members[ordered[i]]):
            j += 1
        run = ordered[i:j]
        run.sort(key=lambda c: (hash((seed, c)) & 0xFFFFFFFF, c))
        out.extend(run)
        i = j
    return out
