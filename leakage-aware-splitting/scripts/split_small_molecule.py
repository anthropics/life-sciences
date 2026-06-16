#!/usr/bin/env python3
"""Small-molecule clustering: Bemis-Murcko scaffolds, ECFP4 fingerprints, Butina.

Scaffold splitting (assign whole Bemis-Murcko scaffold groups to one split) is the
MoleculeNet/DeepChem baseline and is better than random, but it is still
optimistic: benzene and pyridine count as different scaffolds despite near-identical
chemistry, so similar molecules can leak across splits. Always pair a scaffold split
with an ECFP4 Tanimoto nearest-neighbour check, and prefer Butina/fingerprint
clustering when you want a harder, more realistic split.

Requires RDKit. Functions that need it raise a clear error if it is missing so the
orchestrator can report the dependency rather than crash obscurely.
"""

from __future__ import annotations

from typing import Sequence

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.ML.Cluster import Butina
    _RDKIT = True
except Exception:  # pragma: no cover - exercised only when RDKit absent
    _RDKIT = False

PINNED_RDKIT_NOTE = "validated with RDKit 2024+; pin your version in the manifest"


def rdkit_available() -> bool:
    return _RDKIT


def _require_rdkit():
    if not _RDKIT:
        raise RuntimeError("RDKit not installed; `pip install rdkit` for small-molecule splitting")


def bemis_murcko_scaffold(smiles: str) -> str:
    """Canonical Bemis-Murcko scaffold SMILES.

    Returns "" for molecules with no ring system (acyclic) - RDKit's documented
    behaviour - and for invalid SMILES. Callers must treat the empty scaffold as its
    own group, not silently merge all acyclic molecules together (that would be a
    leakage of its own).
    """
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return ""


def scaffold_labels(smiles_list: Sequence[str]) -> list[int]:
    """Cluster labels by Bemis-Murcko scaffold.

    Each distinct scaffold is a cluster. Acyclic molecules (empty scaffold) and
    invalid SMILES each get their OWN singleton cluster rather than being pooled, so
    they are never falsely treated as similar.
    """
    _require_rdkit()
    label_of_scaffold: dict[str, int] = {}
    labels = []
    next_singleton = -1
    for smi in smiles_list:
        scaf = bemis_murcko_scaffold(smi)
        if scaf == "":
            # unique singleton cluster per empty/invalid molecule
            labels.append(next_singleton)
            next_singleton -= 1
            continue
        if scaf not in label_of_scaffold:
            label_of_scaffold[scaf] = len(label_of_scaffold)
        labels.append(label_of_scaffold[scaf])
    # Renumber to contiguous non-negative labels deterministically.
    remap: dict[int, int] = {}
    out = []
    for lbl in labels:
        if lbl not in remap:
            remap[lbl] = len(remap)
        out.append(remap[lbl])
    return out


def morgan_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048):
    """ECFP4-style Morgan fingerprint (radius 2). None for invalid SMILES."""
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def tanimoto(smiles_a: str, smiles_b: str) -> float:
    """ECFP4 Tanimoto similarity in [0, 1]; 0 if either SMILES is invalid."""
    fa = morgan_fingerprint(smiles_a)
    fb = morgan_fingerprint(smiles_b)
    if fa is None or fb is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fa, fb)


def butina_labels(smiles_list: Sequence[str], cutoff: float = 0.4) -> list[int]:
    """Butina clustering on ECFP4 fingerprints; distance cutoff = 1 - Tanimoto.

    cutoff is a *distance*: 0.4 distance == 0.6 Tanimoto similarity merges into a
    cluster. Lower the cutoff (e.g. 0.6 distance / 0.4 similarity) for a harder
    split that separates more near-neighbours. O(n^2); fine for benchmark-scale data.
    """
    _require_rdkit()
    fps = [morgan_fingerprint(s) for s in smiles_list]
    n = len(fps)
    # Condensed distance matrix (lower triangle) as Butina expects.
    dists = []
    for i in range(1, n):
        for j in range(i):
            if fps[i] is None or fps[j] is None:
                dists.append(1.0)
            else:
                dists.append(1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    labels = [0] * n
    for cid, members in enumerate(clusters):
        for m in members:
            labels[m] = cid
    return labels
