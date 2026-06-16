#!/usr/bin/env python3
"""Sequence-modality clustering: produce cluster labels for protein/nucleotide sequences.

Two backends:

* `mmseqs_cluster` - the production path. Shells out to a pinned MMseqs2 binary
  (`easy-cluster`) at a configurable identity/coverage. Default recipe follows the
  methodology survey: `--min-seq-id 0.3 -c 0.8 --cov-mode 0`, plus
  `--cluster-mode 1` (connected component / single linkage) which is what actually
  guarantees between-cluster dissimilarity. Requires MMseqs2 on PATH.

* `kmer_similar_pairs` + clusters_from_pairs - a dependency-free fallback used when
  no binary is available and for hermetic tests. It approximates sequence identity
  with k-mer Jaccard similarity, which is monotonic with identity and good enough
  to demonstrate and test the splitting logic. It is NOT a substitute for MMseqs2
  on real data; the skill warns the user when it falls back.

Both ultimately feed `core.clusters_from_pairs` / cluster labels into the shared
split-assignment core.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Sequence

# The MMseqs2 version the test suite was validated against (see environment.yml).
# Clustering results shift across releases, so the provenance manifest records the
# ACTUAL installed version at run time (see mmseqs_version) and warns when it differs
# from this validated one - never assume the pin is what ran.
PINNED_MMSEQS_VERSION = "18.8cc5c"


def kmer_set(seq: str, k: int = 3) -> frozenset[str]:
    seq = seq.strip().upper()
    if len(seq) < k:
        return frozenset({seq}) if seq else frozenset()
    return frozenset(seq[i:i + k] for i in range(len(seq) - k + 1))


def kmer_jaccard(a: str, b: str, k: int = 3) -> float:
    """Jaccard similarity of k-mer sets; a fast, monotonic proxy for sequence identity."""
    sa, sb = kmer_set(a, k), kmer_set(b, k)
    if not sa and not sb:
        return 1.0
    inter = len(sa & sb)
    union = len(sa | sb)
    return inter / union if union else 0.0


def kmer_similar_pairs(sequences: Sequence[str], threshold: float = 0.3, k: int = 3) -> list[tuple[int, int]]:
    """All pairs with k-mer Jaccard >= threshold. O(n^2); for tests/small data only."""
    pairs = []
    n = len(sequences)
    sets = [kmer_set(s, k) for s in sequences]
    for i in range(n):
        si = sets[i]
        for j in range(i + 1, n):
            sj = sets[j]
            if not si and not sj:
                pairs.append((i, j))
                continue
            union = len(si | sj)
            if union == 0:
                continue
            if len(si & sj) / union >= threshold:
                pairs.append((i, j))
    return pairs


def mmseqs_available() -> bool:
    return shutil.which("mmseqs") is not None


def mmseqs_version() -> str | None:
    """Actual version of the mmseqs binary on PATH (e.g. '18.8cc5c'), or None.

    Recorded in the provenance manifest so the split is auditable against the exact
    binary that produced it, not a hardcoded assumption.
    """
    if not mmseqs_available():
        return None
    try:
        out = subprocess.run(["mmseqs", "version"], capture_output=True, text=True, check=True)
        return out.stdout.strip() or None
    except (subprocess.SubprocessError, OSError):
        return None


def mmseqs_cluster(
    ids: Sequence[str],
    sequences: Sequence[str],
    min_seq_id: float = 0.3,
    coverage: float = 0.8,
    cov_mode: int = 0,
    cluster_mode: int = 1,
    tmp_dir: str | None = None,
) -> list[int]:
    """Cluster sequences with MMseqs2 easy-cluster and return per-item cluster labels.

    cluster_mode defaults to 1 (connected component / single linkage) for the
    strongest leakage guarantee; pass 0 (greedy set cover, the MMseqs2 default) only
    if you understand it does not guarantee between-cluster dissimilarity.

    Raises RuntimeError if the binary is missing so the caller can fall back to the
    k-mer path with an explicit warning rather than silently degrading.
    """
    if not mmseqs_available():
        raise RuntimeError("mmseqs not found on PATH; install MMseqs2 or use the k-mer fallback")

    work = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp(prefix="bioevalsplit_mmseqs_"))
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / "input.fasta"
    with fasta.open("w") as fh:
        for sid, seq in zip(ids, sequences):
            fh.write(f">{sid}\n{seq}\n")

    out_prefix = work / "clu"
    cmd = [
        "mmseqs", "easy-cluster", str(fasta), str(out_prefix), str(work / "tmp"),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--cluster-mode", str(cluster_mode),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    # easy-cluster writes <prefix>_cluster.tsv: representative<TAB>member, one per line.
    rep_of: dict[str, str] = {}
    with (Path(str(out_prefix) + "_cluster.tsv")).open() as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            rep_of[member] = rep

    rep_to_label: dict[str, int] = {}
    labels = []
    for sid in ids:
        rep = rep_of.get(sid, sid)  # singletons may be their own rep
        if rep not in rep_to_label:
            rep_to_label[rep] = len(rep_to_label)
        labels.append(rep_to_label[rep])
    return labels
