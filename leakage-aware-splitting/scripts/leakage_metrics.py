#!/usr/bin/env python3
"""Leakage metrics and split invariants, shared across all modalities.

Two kinds of checks live here:

1. Deterministic invariants any correct split MUST satisfy regardless of data
   (partition completeness/disjointness, cluster integrity, threshold compliance,
   temporal monotonicity, reproducibility). These return structured results so
   tests and the orchestrator can assert on them.

2. Quality statistics that quantify how much leakage remains (the realized max
   train<->test similarity, and the nearest-neighbour similarity distribution).
   These take a `similarity` callable so the same code serves sequences
   (identity), structures (TM-score), and molecules (Tanimoto).

A similarity function here returns a value in [0, 1] where 1 == identical and the
leakage threshold is an upper bound the test set must stay below relative to train.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

SimilarityFn = Callable[[object, object], float]


@dataclass
class InvariantReport:
    name: str
    passed: bool
    detail: str = ""

    def __bool__(self) -> bool:
        return self.passed


# --------------------------------------------------------------------------- #
# Deterministic invariants                                                     #
# --------------------------------------------------------------------------- #

def check_partition(assignment: Sequence[str], n_items: int) -> InvariantReport:
    """Every item is in exactly one split and the union covers all items."""
    if len(assignment) != n_items:
        return InvariantReport("partition_completeness", False,
                               f"assignment length {len(assignment)} != n_items {n_items}")
    missing = [i for i, s in enumerate(assignment) if s is None]
    if missing:
        return InvariantReport("partition_completeness", False,
                               f"{len(missing)} items unassigned (e.g. index {missing[0]})")
    return InvariantReport("partition_completeness", True,
                           f"all {n_items} items assigned to exactly one split")


def check_cluster_integrity(assignment: Sequence[str], cluster_labels: Sequence[int]) -> InvariantReport:
    """No cluster id appears in more than one split (group-wise leakage)."""
    cluster_to_splits: dict[int, set[str]] = {}
    for s, c in zip(assignment, cluster_labels):
        cluster_to_splits.setdefault(c, set()).add(s)
    spanning = {c: sp for c, sp in cluster_to_splits.items() if len(sp) > 1}
    if spanning:
        c0 = next(iter(spanning))
        return InvariantReport("cluster_integrity", False,
                               f"{len(spanning)} clusters span >1 split (e.g. cluster {c0} -> {sorted(spanning[c0])})")
    return InvariantReport("cluster_integrity", True,
                           f"all {len(cluster_to_splits)} clusters confined to one split each")


def check_group_integrity(assignment: Sequence[str], group_keys: Sequence) -> InvariantReport:
    """No non-null group key (patient/site/batch/...) appears in more than one split.

    Same shape as ``check_cluster_integrity`` but on the *raw user-facing keys*, so a
    failure names the leaking patient/site rather than an internal cluster id. Null
    keys are ungrouped singletons and exempt. Passes by construction whenever the
    split was built by composing group edges into the union-find (groups only ever
    merge, never split); checking and reporting it anyway is the diagnostic-instrument
    discipline - the instrument states what it guarantees instead of assuming it.
    """
    key_to_splits: dict[object, set[str]] = {}
    for s, k in zip(assignment, group_keys):
        if k is None:
            continue
        key_to_splits.setdefault(k, set()).add(s)
    spanning = {k: sp for k, sp in key_to_splits.items() if len(sp) > 1}
    if spanning:
        k0 = sorted(spanning, key=lambda k: str(k))[0]
        return InvariantReport("group_integrity", False,
                               f"{len(spanning)} group key(s) span >1 split "
                               f"(e.g. '{k0}' -> {sorted(spanning[k0])})")
    return InvariantReport("group_integrity", True,
                           f"all {len(key_to_splits)} group keys confined to one split each")


def check_temporal_monotonicity(assignment: Sequence[str], dates: Sequence, cutoff_order=("train", "val", "test")) -> InvariantReport:
    """max(date) of an earlier split <= min(date) of a later split.

    `dates` may be anything orderable (ISO strings, datetimes, ints). Splits not
    present in `assignment` are skipped.
    """
    by_split: dict[str, list] = {}
    for s, d in zip(assignment, dates):
        by_split.setdefault(s, []).append(d)
    present = [s for s in cutoff_order if s in by_split]
    for earlier, later in zip(present, present[1:]):
        max_e = max(by_split[earlier])
        min_l = min(by_split[later])
        if max_e > min_l:
            return InvariantReport("temporal_monotonicity", False,
                                   f"max({earlier})={max_e} > min({later})={min_l}")
    return InvariantReport("temporal_monotonicity", True,
                           f"temporal order respected across {present}")


def check_ratio(realized: dict[str, float], requested: dict[str, float], tol: float = 0.15) -> InvariantReport:
    """Realized item-level fractions are within `tol` of requested.

    Tolerance is generous by default because cluster granularity makes exact
    fractions impossible; the point is to catch gross failures, not to demand
    impossible precision. Callers wanting tighter bounds pass a smaller tol.
    """
    worst = 0.0
    worst_split = None
    for s, want in requested.items():
        got = realized.get(s, 0.0)
        if abs(got - want) > worst:
            worst, worst_split = abs(got - want), s
    ok = worst <= tol
    return InvariantReport("ratio_within_tolerance", ok,
                           f"worst deviation {worst:.3f} on '{worst_split}' (tol {tol})")


def check_no_empty_splits(assignment: Sequence[str], requested_fractions: dict) -> InvariantReport:
    """Every split requested with a non-zero fraction actually received items.

    An empty test (or val) set is the sharpest false-reassurance failure mode: any
    metric computed on it is vacuous, yet a generous ratio tolerance can let it pass
    unnoticed. So this is a hard invariant, not a soft warning. With whole-cluster
    assignment it can only happen when there are fewer clusters/groups than requested
    splits (an unsatisfiable request); the fix is fewer splits, a different ratio, or
    more independent units - never silently rebalancing (that would break cluster
    integrity).
    """
    counts = {s: 0 for s in requested_fractions}
    for s in assignment:
        if s in counts:
            counts[s] += 1
    empty = [s for s, f in requested_fractions.items() if f > 0 and counts.get(s, 0) == 0]
    if empty:
        return InvariantReport("no_empty_splits", False,
                               f"requested split(s) {empty} are EMPTY - too few clusters/groups "
                               "to fill all splits at this ratio; reduce the number of splits, "
                               "change the ratio, or gather more independent units")
    return InvariantReport("no_empty_splits", True, "every requested split received items")


def check_reproducible(assignment_a: Sequence[str], assignment_b: Sequence[str]) -> InvariantReport:
    """Two runs with the same seed/input produce identical assignments."""
    ok = list(assignment_a) == list(assignment_b)
    return InvariantReport("reproducibility", ok,
                           "identical" if ok else "assignments differ across identical-seed runs")


# --------------------------------------------------------------------------- #
# Quality statistics (similarity-based)                                        #
# --------------------------------------------------------------------------- #

@dataclass
class LeakageStats:
    max_cross_similarity: float
    threshold: float
    n_violations: int
    nn_similarity: list[float]  # per test item: max similarity to any train item

    @property
    def passes_threshold(self) -> bool:
        return self.max_cross_similarity < self.threshold


def nearest_neighbour_similarity(
    test_items: Sequence,
    train_items: Sequence,
    similarity: SimilarityFn,
) -> list[float]:
    """For each test item, the max similarity to any train item.

    O(len(test) * len(train)); fine for golden/test datasets. Production wrappers
    use MMseqs2/Foldseek search instead of this all-vs-all loop for large sets.
    """
    out = []
    for t in test_items:
        best = 0.0
        for tr in train_items:
            s = similarity(t, tr)
            if s > best:
                best = s
        out.append(best)
    return out


def compute_leakage_stats(
    items: Sequence,
    assignment: Sequence[str],
    similarity: SimilarityFn,
    threshold: float,
    test_split: str = "test",
    train_split: str = "train",
) -> LeakageStats:
    """Realized cross-split leakage: how close is the nearest train item to each test item."""
    test_items = [it for it, s in zip(items, assignment) if s == test_split]
    train_items = [it for it, s in zip(items, assignment) if s == train_split]
    nn = nearest_neighbour_similarity(test_items, train_items, similarity)
    max_cross = max(nn) if nn else 0.0
    n_viol = sum(1 for v in nn if v >= threshold)
    return LeakageStats(max_cross, threshold, n_viol, nn)


def check_scaffold_disjoint(scaffolds: Sequence[str], assignment: Sequence[str]) -> InvariantReport:
    """Bemis-Murcko scaffold sets are disjoint across splits.

    Empty-string scaffolds (acyclic/invalid molecules) are ignored here because each is
    treated as its own singleton upstream; pooling them would create a false overlap.
    """
    by_split: dict[str, set[str]] = {}
    for scaf, s in zip(scaffolds, assignment):
        if scaf == "":
            continue
        by_split.setdefault(s, set()).add(scaf)
    splits = list(by_split)
    for i in range(len(splits)):
        for j in range(i + 1, len(splits)):
            overlap = by_split[splits[i]] & by_split[splits[j]]
            if overlap:
                return InvariantReport("scaffold_disjoint", False,
                                       f"{len(overlap)} scaffolds shared by {splits[i]} & {splits[j]}")
    return InvariantReport("scaffold_disjoint", True, "scaffold sets disjoint across splits")


def ave_bias(
    actives: Sequence,
    inactives: Sequence,
    active_assignment: Sequence[str],
    inactive_assignment: Sequence[str],
    similarity: SimilarityFn,
    test_split: str = "test",
    train_split: str = "train",
) -> float:
    """Asymmetric Validation Embedding (AVE) bias for actives/inactives.

    Quantifies train-validation redundancy that lets a model win by memorizing
    nearest neighbours rather than learning (Wallach & Heifets 2018). Defined as

        AVE = (AA - AI) + (II - IA)

    where, for each test compound, we take the fraction-style nearest-neighbour
    closeness to the train actives (A*) and train inactives (I*), averaged. Values
    near 0 mean a well-balanced (low-bias) split; large positive values mean test
    compounds sit suspiciously close to same-class train compounds. A lightweight
    nearest-neighbour formulation; see the AVE paper for the full embedding version.
    """
    def mean_nn(test_items, train_items):
        if not test_items or not train_items:
            return 0.0
        vals = nearest_neighbour_similarity(test_items, train_items, similarity)
        return sum(vals) / len(vals)

    test_act = [a for a, s in zip(actives, active_assignment) if s == test_split]
    train_act = [a for a, s in zip(actives, active_assignment) if s == train_split]
    test_ina = [a for a, s in zip(inactives, inactive_assignment) if s == test_split]
    train_ina = [a for a, s in zip(inactives, inactive_assignment) if s == train_split]

    aa = mean_nn(test_act, train_act)
    ai = mean_nn(test_act, train_ina)
    ii = mean_nn(test_ina, train_ina)
    ia = mean_nn(test_ina, train_act)
    return (aa - ai) + (ii - ia)
