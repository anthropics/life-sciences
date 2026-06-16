#!/usr/bin/env python3
"""Calibration experiments: prove the split READS generalization honestly.

The split is a diagnostic instrument. These experiments calibrate it against
known truth, in both directions:

* generalization_gap: train ONE cheap fixed model on a naive random split and on
  the leakage-aware split. When leakage is real, the random split conceals a gap
  (inflated score) that the leakage-aware split exposes (realistic, lower score).
  When the signal is genuinely generalizable, BOTH splits score well and there is
  no gap - so a large gap is evidence of leakage, not of a broken model.

* nn_baseline_collapse: the k-NN "memorize the nearest training example" baseline
  scores well on a leaky split and collapses toward the mean-label baseline on an
  honest split.

Uses a k-mer featurizer for sequences and sklearn RandomForest as the fixed model,
so it is cheap, deterministic, and dependency-light.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from sklearn.ensemble import RandomForestClassifier

from .split import make_split


def kmer_features(sequences: Sequence[str], k: int = 2) -> np.ndarray:
    """Fixed k-mer count feature matrix. k=2 over the 20 amino acids -> 400 features."""
    AA = "ACDEFGHIKLMNPQRSTVWY"
    kmers = [a + b for a in AA for b in AA] if k == 2 else None
    if kmers is None:
        raise ValueError("kmer_features currently supports k=2")
    index = {km: i for i, km in enumerate(kmers)}
    X = np.zeros((len(sequences), len(kmers)), dtype=float)
    for r, seq in enumerate(sequences):
        seq = seq.upper()
        for i in range(len(seq) - 1):
            j = index.get(seq[i:i + 2])
            if j is not None:
                X[r, j] += 1.0
    return X


def _random_assignment(n: int, seed: int, frac_test: float = 0.25) -> list[str]:
    rng = random.Random(seed)
    return ["test" if rng.random() < frac_test else "train" for _ in range(n)]


def _accuracy(model, X_test, y_test) -> float:
    if len(y_test) == 0:
        return float("nan")
    return float((model.predict(X_test) == np.asarray(y_test)).mean())


@dataclass
class GapResult:
    random_split_score: float
    leakage_aware_score: float
    gap: float  # random - leakage_aware; large positive == leakage was concealing weakness


def generalization_gap(
    ids: Sequence[str],
    sequences: Sequence[str],
    labels: Sequence[int],
    threshold: float = 0.3,
    seed: int = 0,
) -> GapResult:
    """Compare a fixed model's test accuracy under a random vs leakage-aware split."""
    X = kmer_features(sequences)
    y = np.asarray(labels)

    # --- random split (the naive, leakage-prone baseline) ---
    rand = _random_assignment(len(sequences), seed)
    tr = [i for i, s in enumerate(rand) if s == "train"]
    te = [i for i, s in enumerate(rand) if s == "test"]
    m1 = RandomForestClassifier(n_estimators=100, random_state=seed)
    m1.fit(X[tr], y[tr])
    rand_score = _accuracy(m1, X[te], y[te])

    # --- leakage-aware split from the skill ---
    report = make_split(ids, sequences, modality="sequence",
                        fractions={"train": 0.75, "test": 0.25},
                        threshold=threshold, seed=seed, prefer_tool=False)
    a = report["assignment"]
    tr2 = [i for i, s in enumerate(a) if s == "train"]
    te2 = [i for i, s in enumerate(a) if s == "test"]
    m2 = RandomForestClassifier(n_estimators=100, random_state=seed)
    m2.fit(X[tr2], y[tr2])
    aware_score = _accuracy(m2, X[te2], y[te2])

    return GapResult(rand_score, aware_score, rand_score - aware_score)
