# Group / metadata-axis splitting (the second leakage axis)

Biological leakage has **two orthogonal axes**. The similarity axis (sequence,
structure, scaffold) is handled by the per-modality references. This one is the
**group axis**: replicate units that share a *source entity* must travel to the same
split, or the model is scored on memorizing the entity instead of learning the task.

The canonical footgun is histology: a whole-slide image is cut into thousands of
tiles, slides come from patients, patients from sites/scanners. Split by **tile** and
the same patient lands in both train and test; the model memorizes that patient's
staining and tissue and reports inflated accuracy. The unit you care about
generalizing to is the **patient**, not the tile. The same shape recurs across
biology:

| Domain | Rows (what you split) | Group (what must stay together) |
|--------|-----------------------|----------------------------------|
| Digital pathology | tiles / patches | slide ⊂ patient ⊂ site/scanner |
| Clinical | samples / time-points | subject / donor / case |
| High-throughput screening, imaging | wells / images | plate / batch / run / site |
| Sequencing | reads / libraries | sample / specimen / replicate |
| Pooled datasets | examples | study / source / dataset / cohort |
| Comparative genomics | genes / proteins | organism / strain |

## The orienting question

> **What is the unit of generalization?** At deployment, what will be *new* that the
> model has never seen - a new patient? a new site? a new study?

Split so that whole unit is unseen in test. That single question picks the grouping
column. When the data nests (tile ⊂ slide ⊂ patient ⊂ site), the **coarsest** unit
gives the most protection - but coarser also means fewer groups, so an 80/10/10 may
become unsatisfiable (you cannot put 3 sites into 3 splits and keep a test set). The
skill surfaces the whole chain and the satisfiability of each level so you can choose.

## Orthogonal to, and composable with, similarity

A group key is just another source of "these items must co-cluster" edges. The skill
emits one spanning edge per group member and feeds them into the *same* single-linkage
union-find that the similarity tools feed (`core.clusters_from_pairs`). Composition is
therefore one pass over the union of both edge sets:

```
similarity edges (MMseqs2/Foldseek/scaffold)  ∪  group edges (patient/site/batch)
        → connected components → whole-component assignment to splits
```

Every existing guarantee (whole-component assignment, cluster integrity,
reproducibility, stratification) carries over unchanged. You can also run the group
axis **alone** (`--modality metadata`) when no biological-similarity tool applies -
e.g. histology, where the rows are images, not sequences.

## Detection is propose-then-confirm (never silent)

**Automatic grouping-column detection is not a solved problem** - no off-the-shelf
tool does it reliably, because whether a column is a generalization unit is a question
about your task, not about the data. So the skill **profiles and proposes; it never
auto-applies.** `scripts/group_detect.py:propose_grouping` is pure: it returns a
ranked proposal and changes nothing. The CLI enforces the discipline - with
`--auto-detect-groups` it prints the ranking and **exits without splitting** unless
you also pass a confirmed `--group-col`.

The proposal ranks each column by: dictionary match (name families like
patient/site/batch, weighted toward true generalization units), hierarchy coarseness
(coarser = more protection), and satisfiability for the requested ratio. It **refuses
to recommend** (disqualifies, with a stated reason) columns that are unique-per-row
(grouping is a no-op - it is the row id), constant (one group cannot be split),
near-constant (one group holds ~all rows), or equal to the stratification/label column
(grouping by the target leaks it). Ties at the top are surfaced explicitly so you ask
rather than guess.

## CLI walkthrough (histology)

```bash
# 1) Profile and propose - prints the ranking, splits nothing.
python -m scripts.split --modality metadata --metadata tiles.csv \
    --id-col tile_id --auto-detect-groups
#   RECOMMENDED: --group-col patient_id  (15 groups ...)
#    - patient_id   score=2.50   satisfiable
#    - slide_id     score=1.40   satisfiable
#    - site         score=...    NOT satisfiable (only 3 groups)
#    - tile_id      DISQUALIFIED (unique per row)

# 2) Confirm the unit of generalization and split.
python -m scripts.split --modality metadata --metadata tiles.csv \
    --id-col tile_id --group-col patient_id --train 0.8 --val 0.1 --test 0.1 --out splits/
```

Compose with a sequence split by adding `--metadata`/`--group-col` to a sequence run:

```bash
python -m scripts.split --modality sequence --fasta seqs.fasta \
    --metadata meta.csv --id-col seq_id --group-col donor_id --out splits/
```

## Hierarchies (pick the coarsest you can afford)

When columns nest, the skill recovers the chain with functional dependencies (A
determines B iff every value of A co-occurs with exactly one value of B) and scores
**coarseness** = how many finer columns determine it. Recommend the coarsest, but
respect satisfiability: with only a handful of sites, splitting by site empties a
split - then patient is the right compromise, and the report says so.

## Stratified group splits

To balance a class label *and* keep groups intact, compose first (group edges →
components) then assign with `stratify_labels` (each component's majority class is its
stratum). Balance is bounded by group granularity: if every member of a class lives in
one or two groups, you cannot finely balance it - that is a property of the data, and
the realized ratio is reported rather than faked.

## Failure modes and refusals

| Symptom | What the skill does |
|---------|---------------------|
| Grouping column is unique per row | Disqualify: grouping is a no-op (it is the row id). |
| One group holds ~all rows (near-constant) | Refuse to recommend: it would force ~all data into one split. |
| Fewer groups than splits | A split is unavoidably empty; the `no_empty_splits` invariant FAILS (hard). The proposal flags `satisfiable_for = False` up front. |
| Enough groups but ratio not achievable (e.g. 4 groups, 80/10/10) | The packer seeds the small splits first so none is empty, but `ratio_within_tolerance` fails - the honest "cannot hit this ratio" signal. |
| Ambiguous (several equally-good keys) | Return all tied candidates; ask the user. |
| Group column equals the label/stratify column | Disqualify: grouping by the target leaks it. |
| Null group values | Treated as ungrouped singletons; counted in `n_null_group_keys`. Use `--drop-null-groups` to exclude them. |
| Requested ratio cannot be met with whole groups | Never silently rebalance (that would break group integrity); report the realized ratio and warn. |

## Interpreting the report

The `grouping` block reports `n_groups`, the group-size distribution, `n_null_group_keys`,
`composed_with_similarity`, and **`clusters_merged_across_groups`** - the number of
final clusters that contain more than one original group key. That last number is the
honest "did biological similarity glue two groups together?" measure: if two patients
share a near-identical sequence, they must co-travel, and your effective number of
independent units is smaller than your patient count. `check_group_integrity` then
confirms (and names any offender) that no raw group key spans two splits.

## Citations

- Group/blocked cross-validation: Roberts et al., "Cross-validation strategies for
  data with temporal, spatial, hierarchical, or phylogenetic structure," Ecography 2017.
- `GroupKFold` / `StratifiedGroupKFold`: Pedregosa et al., scikit-learn, JMLR 2011.
- Patient-level leakage in medical imaging: Oakden-Rayner et al., "Hidden
  stratification," ACM CHIL 2020; Kapoor & Narayanan, "Leakage and the reproducibility
  crisis in ML-based science," Patterns 2023 (data leakage taxonomy, incl. group leakage).
- Batch effects as leakage: Leek et al., "Tackling the widespread and critical impact
  of batch effects," Nat. Rev. Genet. 2010.
