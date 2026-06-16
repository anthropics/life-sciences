#!/usr/bin/env python3
"""Temporal / deposition-date splitting: train on the past, test on the future.

Time-split is the closest cheap proxy to true prospective performance (Sheridan
2013): random selection is too optimistic, leave-class-out too pessimistic. AlphaFold2
validated this way ("all structures deposited after our training cutoff"). It mimics
real deployment and captures all leakage types as they actually emerge.

Caveat encoded here: a temporal split alone does NOT guarantee novelty - the
post-cutoff set can still contain close homologs of pre-cutoff entries. Best practice
combines a temporal cutoff with a similarity filter; this module produces the
temporal assignment and reports which test items still have a near-duplicate in train
so the caller can additionally filter them.
"""

from __future__ import annotations

from typing import Callable, Sequence


def temporal_split(
    dates: Sequence,
    fractions: dict[str, float] | None = None,
) -> list[str]:
    """Assign splits by chronological order so max(train) <= min(val) <= min(test).

    `dates` may be ISO strings, datetimes, or ints (anything orderable). Items are
    sorted by date; the earliest `train` fraction become train, then val, then the
    latest become test. Ties on the boundary date are resolved by keeping all items
    with the same date on the earlier side, which can shift realized fractions - that
    is the honest behaviour (you cannot place identical-date items on both sides
    without leaking the cutoff).
    """
    if fractions is None:
        fractions = {"train": 0.8, "val": 0.1, "test": 0.1}
    order = [s for s in ("train", "val", "test") if fractions.get(s, 0) > 0]
    n = len(dates)
    idx_sorted = sorted(range(n), key=lambda i: dates[i])

    # Cumulative cut points by fraction.
    assignment = [None] * n
    cuts = []
    acc = 0.0
    for s in order[:-1]:
        acc += fractions[s]
        cuts.append(int(round(acc * n)))
    bounds = [0] + cuts + [n]

    for split, lo, hi in zip(order, bounds, bounds[1:]):
        for pos in range(lo, hi):
            assignment[idx_sorted[pos]] = split

    # Enforce no boundary-date straddling: if the item at a boundary shares a date
    # with the previous split's last item, push it back to the earlier split.
    for b in range(1, len(order)):
        boundary_pos = bounds[b]
        if 0 < boundary_pos < n:
            d_before = dates[idx_sorted[boundary_pos - 1]]
            pos = boundary_pos
            while pos < n and dates[idx_sorted[pos]] == d_before:
                assignment[idx_sorted[pos]] = order[b - 1]
                pos += 1
    return assignment


def residual_homology_in_temporal_split(
    items: Sequence,
    assignment: Sequence[str],
    similarity: Callable[[object, object], float],
    threshold: float,
) -> list[int]:
    """Test-item indices that still have a >=threshold neighbour in train.

    A temporal split can pass chronology yet leak by homology. This flags the test
    items the caller should additionally filter (or at least disclose).
    """
    test_idx = [i for i, s in enumerate(assignment) if s == "test"]
    train_items = [items[i] for i, s in enumerate(assignment) if s == "train"]
    flagged = []
    for i in test_idx:
        if any(similarity(items[i], tr) >= threshold for tr in train_items):
            flagged.append(i)
    return flagged
