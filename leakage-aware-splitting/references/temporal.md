# Temporal / deposition-date splitting

Train on data from before a cutoff date, test on data after it. The closest cheap
proxy to true prospective performance (Sheridan 2013): random selection is too
optimistic, leave-class-out too pessimistic. AlphaFold2 validated this way ("all
structures deposited after our training cutoff"); CASP (biennial blind) and CAMEO
(continuous) are the gold-standard temporal benchmarks.

## Strengths

- Mimics real prospective deployment.
- Captures all leakage types simultaneously (sequence, structure, metadata) as they
  actually emerge over time.
- No arbitrary similarity threshold to choose.

## Weaknesses (encoded in the tooling)

- The post-cutoff set may still contain close homologs of pre-cutoff entries (the same
  family keeps getting deposited), so a temporal split alone does **not** guarantee
  novelty. `scripts/split_temporal.py:residual_homology_in_temporal_split` flags test
  items that still have a near-duplicate in train so you can additionally filter or
  disclose them. **Best practice = temporal cutoff + similarity filter.**
- Can be "too easy" vs an identity split, or "too strict" if real deployment is on
  same-family inputs.
- Requires reliable timestamps; the test set shrinks as the cutoff moves.

## Recipe

`temporal_split(dates, fractions)` sorts items chronologically and assigns the earliest
to train, then val, then the latest to test, guaranteeing
max(train) <= min(val) <= min(test). Items sharing the boundary date are kept on the
earlier side (you cannot place identical-date items on both sides without leaking the
cutoff), which can shift realized fractions - that is honest behaviour, reported, not
hidden.

## When no dates exist

**SIMPD** (Landrum et al., J. Cheminform. 2023) uses a multi-objective genetic algorithm
to generate simulated time-splits mimicking real medicinal-chemistry train/test
differences (99 public ChEMBL-derived sets), a useful golden-reference comparator.

## Citations

- Sheridan, time-split validation, JCIM 2013.
- Jumper et al., AlphaFold2, Nature 2021 (post-cutoff validation).
- Landrum et al., SIMPD, J. Cheminform. 2023.
