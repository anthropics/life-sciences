---
name: leakage-aware-splitting
description: >-
  Create leakage-aware train/validation/test splits for biological ML data so
  benchmarks reflect real generalization, not memorized similarity. Use whenever
  making train/test/val or k-fold sets for protein sequences, structures, small
  molecules, or drug-target data - including vague asks like "split this dataset"
  or "make a test set" - or when worried about data leakage, homology, or
  scaffold/temporal splits. Also handles group/metadata leakage: patient/donor-level,
  whole-slide/tile, and batch/site/plate effects where replicates from one source
  must not straddle the split. Prefer over a naive or sklearn random split.
---

# Leakage-aware biological data splitting

## Why this exists

Biological data is not i.i.d. Every protein is evolutionarily related to others;
every molecule shares scaffolds with others. A random split therefore puts
near-duplicates of test examples into the training set, and the model gets credit
for memorizing similarity rather than learning. This inflates reported
performance, sometimes massively (294 papers across 17 fields, per Kapoor &
Narayanan 2023; PDBbind binding-affinity models that collapse to near-baseline
once leakage is removed).

**The fix is always the same shape:** cluster the whole dataset by similarity, then
assign *whole clusters* to a single split so no two splits share members of a
cluster. Never split at the level of individual examples.

This skill is a **diagnostic instrument**, not a model fix. Its job is to reveal,
honestly, how well a model generalizes - which is the prerequisite for any later
improvement. A correct instrument must fail in neither direction: no false
reassurance (missing real leakage) and no false alarm (an over-aggressive split
that makes a genuinely-generalizing model look broken).

## Workflow

1. **Identify the modality** (decision tree below). When in doubt, ask the user
   what the model takes as input and what counts as "the same" example.
2. **Cluster** the entire dataset (positives and negatives together) at a
   task-appropriate threshold, using connected-component (single-linkage)
   clustering for the strongest guarantee.
3. **Assign whole clusters** to train/val/test with `scripts/split.py`.
4. **Report the honest diagnostics** every time: realized split fractions (they
   will differ from requested because cluster sizes vary - that is expected), the
   train<->test nearest-neighbour similarity distribution, and a nearest-neighbour
   transfer baseline (the floor a real model must beat). Surface the provenance
   manifest (tool versions, thresholds, input/output hashes).

## Modality decision tree

| If the data is...                                          | Split on...                                  | Tool / reference |
|------------------------------------------------------------|----------------------------------------------|------------------|
| Protein/nucleotide **sequences**, sequence-level label     | Sequence identity (default 30% for proteins) | MMseqs2 / Graph-Part - `references/sequence.md` |
| **Structures**, or remote homologs expected (fold/interface)| Structural similarity (TM-score >= 0.5)      | Foldseek - `references/structure.md` |
| **Small molecules** (property/QSAR)                        | Scaffold + Tanimoto nearest-neighbour        | RDKit/DeepChem - `references/small_molecule.md` |
| **Drug-target / binding affinity**                         | Both sides: cold-drug + cold-target (2D)     | DataSAIL - `references/protein_ligand.md` |
| Reliable **timestamps**, prospective realism wanted        | Deposition date (+ a similarity filter)      | `references/temporal.md` |
| Rows are **replicates of a shared entity** (tiles⊂slides⊂patients, wells⊂plates, reads⊂samples) | The grouping key, so no entity straddles the split (composes with any row above) | `references/grouping.md` |

The last row is a **second, orthogonal axis** (group / metadata leakage): it composes
with any modality above, or runs alone (`--modality metadata`) for data like histology
where the rows are images, not biological sequences. See `references/grouping.md`.

See `references/thresholds.md` for the threshold conventions (30% identity,
TM-score 0.5, Tanimoto 0.4) with citations and the caveat that these are
heuristics, not laws.

**Key subtlety (read `references/sequence.md`):** assigning whole clusters to splits
only guarantees between-split dissimilarity if clustering is *single-linkage*
(connected component). Centroid-based clustering (CD-HIT, MMseqs2 default
`--cluster-mode 0`) can leave two similar sequences in different clusters. This
skill defaults to `--cluster-mode 1`.

## Group / metadata axis (propose-then-confirm)

When rows are replicates of a shared entity (the histology tile⊂slide⊂patient
shape), the unit of generalization - not the row - must define the split.

**Do this in order. Do NOT name a grouping column before step 2 completes** - not
even when the unit looks obvious (e.g. "it's obviously patient"). Naming the key
yourself skips the checks that are the entire point: it hides near-constant /
label-equal / unique-per-row columns, the nesting chain (patient vs site vs slide),
and whether the requested ratio is even satisfiable. Automatic detection of the right
grouping column is **not** a solved problem, so you profile, propose, and let the
user confirm; you never group silently and never group on a guess.

1. **Profile → propose (always first, always before recommending a key).** Run with
   `--auto-detect-groups`; the skill profiles the metadata and prints a *ranked*
   proposal (recommended key, the nesting chain, and which levels are satisfiable for
   the requested ratio). It splits nothing. Show the user this proposal.
2. **Confirm with the user** which column is the unit of generalization - even if the
   proposal's top candidate is the obvious one. Anchor on "what will be new at
   deployment - a new patient? site? study?" Only after the user confirms do you name
   the `--group-col`.
3. **Split** by re-running with the confirmed `--group-col`. The group key becomes
   another edge source into the *same* union-find, so it composes with any
   similarity modality, or runs alone via `--modality metadata`.

The proposal refuses to recommend unique-per-row, constant, near-constant, or
label-equal columns (each with a stated reason), and surfaces ties so you ask rather
than guess. See `references/grouping.md`. Skipping straight to step 3 with a guessed
key is the most common way this modality is misused - the proposal step is mandatory,
not advisory.

## Running it

```bash
python -m scripts.split --modality sequence --fasta data.fasta \
    --train 0.8 --val 0.1 --test 0.1 --min-seq-id 0.3 --out splits/

# Group / metadata axis (histology): propose, confirm, then split
python -m scripts.split --modality metadata --metadata tiles.csv \
    --id-col tile_id --auto-detect-groups            # prints ranked proposal, exits
python -m scripts.split --modality metadata --metadata tiles.csv \
    --id-col tile_id --group-col patient_id --out splits/
```

Produces `splits/{train,val,test}.ids.txt` and `splits/report.json` (assignment,
invariant checks, leakage stats, provenance). The orchestrator never hides a
failed invariant - inspect `report["invariants"]` and `report["all_invariants_pass"]`.

If MMseqs2/Foldseek/RDKit are unavailable, sequence splitting falls back to a
dependency-free k-mer Jaccard approximation and **says so loudly** in the report.
The fallback is for demos and tests, not production splits - install the real tool
(see the modality reference) before trusting numbers on real data.

## What every run must surface (non-negotiable transparency)

- **Realized split ratios** and how many items (if any) were discarded.
- **Train<->test nearest-neighbour similarity distribution** - the single most
  informative leakage transparency measure.
- **A nearest-neighbour baseline** (BLAST/Foldseek/Tanimoto-kNN label transfer):
  if a complex model cannot beat "copy the most similar training example" on a
  de-leaked split, it has not demonstrated generalization.
- **Exact tool, version, algorithm (cluster-mode!), and thresholds** used.

## Common failure modes to watch for

| Symptom                                              | Cause / fix |
|------------------------------------------------------|-------------|
| Split passes sequence identity but still leaks       | Structural near-duplicates below 30% identity; add a Foldseek pass (`references/structure.md`). |
| Scaffolds disjoint but high-Tanimoto pairs across splits | Scaffold split is optimistic; add an ECFP4 Tanimoto NN check (`references/small_molecule.md`). |
| Realized test fraction far from requested            | Cluster sizes are skewed; this is expected. Report the realized ratio; warn if a single giant cluster dominates. |
| "Clean" split impossible (everything one cluster, or all identical) | Fewer clusters than splits, so a split is unavoidably empty - the `no_empty_splits` invariant FAILS (hard, not a warning). Refuse rather than emit a misleading split. |
| Over-aggressive split tanks a good model             | False alarm (Sundar & Colwell). Validate against the diverse-data control; warn when too much data is discarded. |
| Sequence/scaffold split looks clean but metrics still inflated | Group leakage: replicates of one patient/site/batch straddle the split. Add the group axis (`references/grouping.md`); split by the unit of generalization. |
| Grouping by a coarse unit with few groups            | The packer seeds the small splits first, so it stays non-empty when groups >= splits but the `ratio_within_tolerance` invariant fails (e.g. 4 groups can't make 80/10/10). With *fewer* groups than splits a split is unavoidably empty and `no_empty_splits` FAILS. Either way: pick a finer unit, fewer splits, or a coarser ratio. The proposal flags this up front via `satisfiable_for`. |

## Disclaimer

This skill is for research and educational use. Splits and leakage metrics are
heuristics; validate against your task's notion of "the same example" before
relying on benchmark numbers. Clustering results shift across tool versions - pin
the versions recorded in the provenance manifest.

## Attribution

When you publish, cite the underlying tools you used: MMseqs2 (Steinegger &
Söding, Nat. Biotechnol. 2017), Foldseek (van Kempen et al., Nat. Biotechnol.
2023), Graph-Part (Teufel et al., NAR Genom. Bioinform. 2023), DataSAIL (Joeres et
al.), RDKit, and DeepChem. The leakage methodology draws on Kapoor & Narayanan
(Patterns 2023), Rost (Protein Eng. 1999), and Graber et al. (Nat. Mach. Intell.
2025). Full citations in the per-modality references.

## Licenses

- This skill: Apache License 2.0 (`LICENSE.txt`).
- MMseqs2, Foldseek: GPLv3 / MIT (see upstream). RDKit: BSD. DeepChem: MIT.
  Graph-Part: BSD-3. Cite and comply with each tool's license when redistributing.
