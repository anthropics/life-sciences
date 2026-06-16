# leakage-aware-splitting

A Claude skill for creating **leakage-aware train/validation/test splits** for
biological machine-learning data, so model benchmarks reflect real generalization
instead of memorized similarity.

Biological data is not i.i.d.: every protein is evolutionarily related to others and
every molecule shares scaffolds with others, so a random split puts near-duplicates of
test examples into training and the model gets credit for memorizing similarity. This
skill clusters the whole dataset by similarity and assigns **whole clusters** to a
single split, then reports the honest diagnostics every benchmark should include.

## What's here

| Path | What it is |
|------|-----------|
| `SKILL.md` | The skill: modality decision tree, workflow, transparency requirements. |
| `scripts/` | The implementation. `split.py` orchestrates; `core.py` is the tool-agnostic split-assignment core; per-modality clustering + shared `leakage_metrics.py`, `baselines.py`, `calibration.py`. |
| `references/` | One doc per modality (`sequence`, `structure`, `small_molecule`, `protein_ligand`, `temporal`) plus `grouping.md` and `thresholds.md`. |

## Quick start

Run these from this directory (`leakage-aware-splitting/`):

```bash
# Option A: pip (pure-Python; sequence falls back to a k-mer approximation,
# structure modality unavailable):
pip install -r requirements.txt

# Option B: conda (recommended) - adds the MMseqs2 and Foldseek binaries so the
# real sequence/structure paths run:
conda env create -f environment.yml
conda activate leakage-aware-splitting

# Sequence split (uses MMseqs2 if present, else a k-mer approximation):
python -m scripts.split --modality sequence --fasta data.fasta \
    --train 0.8 --val 0.1 --test 0.1 --min-seq-id 0.3 --out splits/
```

## Modalities

sequence (MMseqs2) · structure (Foldseek) · small_molecule (RDKit scaffold/Butina/
Tanimoto) · protein_ligand (double-cold / DataSAIL) · temporal (deposition date) ·
metadata/grouped (group detection).

## Design principle

The skill is a **diagnostic instrument**, not a model fix. It must fail in neither
direction: no false reassurance (missing real leakage) and no false alarm (an
over-aggressive split that makes a genuinely-generalizing model look broken). The
calibration tests verify both.

License: Apache-2.0. See `SKILL.md` for attribution and per-tool licenses.
