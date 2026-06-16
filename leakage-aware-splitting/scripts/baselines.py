#!/usr/bin/env python3
"""Nearest-neighbour and mean-label baselines: the honest floor a model must beat.

The point of a leakage-aware split is revealed by these baselines. On a leaky split,
"predict the label of the most similar training example" scores competitively
(memorization works). On an honest split it collapses toward the trivial mean-label
baseline. A complex model that cannot beat the kNN baseline on a de-leaked split has
not demonstrated generalization (this is the CAFA BLAST-baseline logic, and the AVE
nearest-neighbour logic for chemistry).

These operate on any modality via a `similarity` callable, so the same code scores
sequences (identity), structures (TM-score), and molecules (Tanimoto).
"""

from __future__ import annotations

from typing import Callable, Sequence

SimilarityFn = Callable[[object, object], float]


def knn_label_transfer(
    test_items: Sequence,
    train_items: Sequence,
    train_labels: Sequence[float],
    similarity: SimilarityFn,
    k: int = 5,
) -> list[float]:
    """Predict each test label as the mean label of its k most-similar train items.

    Works for regression (mean of neighbour values) and binary classification
    (mean -> probability). Returns one prediction per test item. Ties are broken by
    train index for determinism.
    """
    preds = []
    for t in test_items:
        sims = sorted(
            range(len(train_items)),
            key=lambda j: (similarity(t, train_items[j]), -j),
            reverse=True,
        )[:k]
        if not sims:
            preds.append(float("nan"))
            continue
        preds.append(sum(train_labels[j] for j in sims) / len(sims))
    return preds


def mean_label_baseline(train_labels: Sequence[float], n_test: int) -> list[float]:
    """Assign the mean training label to every test item (the trivial floor)."""
    mean = sum(train_labels) / len(train_labels) if train_labels else 0.0
    return [mean] * n_test


def rmse(pred: Sequence[float], true: Sequence[float]) -> float:
    n = len(true)
    if n == 0:
        return float("nan")
    return (sum((p - y) ** 2 for p, y in zip(pred, true)) / n) ** 0.5


def baseline_report(
    items: Sequence,
    labels: Sequence[float],
    assignment: Sequence[str],
    similarity: SimilarityFn,
    k: int = 5,
) -> dict:
    """Compare kNN vs mean-label baseline on the given split.

    Returns RMSEs and the gap. A small gap (kNN barely better than mean-label) on a
    de-leaked split is the healthy signal: similarity memorization buys little, so
    the split is honest. A large gap suggests residual leakage the kNN exploits.
    """
    test_idx = [i for i, s in enumerate(assignment) if s == "test"]
    train_idx = [i for i, s in enumerate(assignment) if s == "train"]
    test_items = [items[i] for i in test_idx]
    test_true = [labels[i] for i in test_idx]
    train_items = [items[i] for i in train_idx]
    train_labels = [labels[i] for i in train_idx]

    knn_pred = knn_label_transfer(test_items, train_items, train_labels, similarity, k=k)
    mean_pred = mean_label_baseline(train_labels, len(test_idx))

    knn_rmse = rmse(knn_pred, test_true)
    mean_rmse = rmse(mean_pred, test_true)
    return {
        "knn_rmse": knn_rmse,
        "mean_label_rmse": mean_rmse,
        "gap": mean_rmse - knn_rmse,  # how much memorization buys; small == honest
        "k": k,
        "n_test": len(test_idx),
    }
