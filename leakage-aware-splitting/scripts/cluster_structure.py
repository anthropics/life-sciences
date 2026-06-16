#!/usr/bin/env python3
"""Structure-modality clustering with Foldseek.

When remote homologs share a fold but not detectable sequence identity, a sequence
split leaks. Foldseek converts structures to a 3Di structural alphabet and searches
via MMseqs2, making all-vs-all structural comparison tractable (TM-align would take
millennia at that scale). Cluster structures, then assign whole structural clusters
to one split. Standard fold-equivalence threshold: TM-score >= 0.5; Foldseek also
uses E-value < 0.01.

Requires the Foldseek binary on PATH. There is no pure-Python fallback for real
structural similarity (it needs 3D coordinates), so when Foldseek is absent the
caller must either supply precomputed cluster labels or fall back to a sequence
split with an explicit warning that structural leakage is unchecked.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Sequence

PINNED_FOLDSEEK_NOTE = "pin your Foldseek version in the manifest; clustering shifts across releases"
# Foldseek version the test suite was validated against (see environment.yml).
PINNED_FOLDSEEK_VERSION = "10.941cd33"


def foldseek_available() -> bool:
    return shutil.which("foldseek") is not None


def foldseek_version() -> str | None:
    """Actual version of the foldseek binary on PATH, or None. Recorded in the
    provenance manifest so the structural split is auditable against the exact binary."""
    if not foldseek_available():
        return None
    try:
        out = subprocess.run(["foldseek", "version"], capture_output=True, text=True, check=True)
        return out.stdout.strip() or None
    except (subprocess.SubprocessError, OSError):
        return None


def foldseek_cluster(
    pdb_paths: Sequence[str],
    ids: Sequence[str] | None = None,
    min_seq_id: float = 0.0,
    coverage: float = 0.8,
    tmscore_threshold: float = 0.5,
    tmp_dir: str | None = None,
) -> list[int]:
    """Cluster structures with `foldseek easy-cluster`; return per-item cluster labels.

    Parameters mirror the structural-leakage recipe: cluster at TM-score >= 0.5 with
    80% coverage. `pdb_paths` are paths to .pdb/.cif files in item order; `ids`
    default to the file stems.

    Raises RuntimeError if Foldseek is missing so the caller can choose an explicit
    fallback rather than silently skipping structural leakage control.
    """
    if not foldseek_available():
        raise RuntimeError("foldseek not found on PATH; install Foldseek for structural splits")
    if ids is None:
        ids = [Path(p).stem for p in pdb_paths]

    work = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp(prefix="bioevalsplit_foldseek_"))
    indir = work / "structures"
    indir.mkdir(parents=True, exist_ok=True)
    # Foldseek easy-cluster takes a directory of structures. Symlink/copy with stable names.
    name_to_id = {}
    for p, sid in zip(pdb_paths, ids):
        src = Path(p)
        dst = indir / src.name
        if not dst.exists():
            dst.write_bytes(src.read_bytes())
        name_to_id[src.stem] = sid

    out_prefix = work / "clu"
    cmd = [
        "foldseek", "easy-cluster", str(indir), str(out_prefix), str(work / "tmp"),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--tmscore-threshold", str(tmscore_threshold),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    rep_of: dict[str, str] = {}
    with (Path(str(out_prefix) + "_cluster.tsv")).open() as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            rep_of[member] = rep

    rep_to_label: dict[str, int] = {}
    labels = []
    for sid_stem in [Path(p).stem for p in pdb_paths]:
        rep = rep_of.get(sid_stem, sid_stem)
        if rep not in rep_to_label:
            rep_to_label[rep] = len(rep_to_label)
        labels.append(rep_to_label[rep])
    return labels
