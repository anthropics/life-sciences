# Sequence-identity splitting (proteins and nucleotides)

The primary, most common case: the model takes a sequence and predicts a
sequence-level property, and there are no structures (or structure is irrelevant
to the task).

## When to use this

- Protein or nucleotide sequences with a sequence-level label (function, signal
  peptide, localization, stability, expression, fitness).
- No structures available, or the task does not depend on fold/interface.
- If remote homologs that share a fold but not detectable sequence identity matter
  for your task, sequence identity is *necessary but not sufficient* - add a
  structural pass (`structure.md`).

## The recipe

1. Cluster the **entire** dataset (positives + negatives together) by sequence
   identity.
2. Assign **whole clusters** to one split each.
3. Verify the realized max train<->test identity is below the configured cutoff and
   report the nearest-neighbour identity distribution.

### Default MMseqs2 command

```bash
mmseqs easy-cluster input.fasta clu tmp \
    --min-seq-id 0.3 -c 0.8 --cov-mode 0 --cluster-mode 1
```

- `--min-seq-id 0.3` - 30% identity, the canonical protein threshold (Rost's
  twilight zone; see `thresholds.md`). Use 0.25 for stricter homology removal
  (TAPE), lower (<0.2) for structure-prediction-grade rigor. Nucleotides need a
  **higher** threshold (random identity baseline ~25% for a 4-letter alphabet).
- `-c 0.8 --cov-mode 0` - bidirectional 80% coverage (divides by the longer
  sequence), so short fragments are not absorbed into clusters of long sequences.
- `--cluster-mode 1` - **connected component / single linkage.** This is the
  important one (see below).

The skill's `scripts/cluster_sequence.py:mmseqs_cluster` runs exactly this and
returns per-item cluster labels for the shared split core.

## The cluster-mode subtlety (do not skip)

Assigning whole clusters to splits guarantees between-split dissimilarity **only if
clustering is single-linkage**. MMseqs2 offers:

- `--cluster-mode 0` greedy set cover (the MMseqs2 default) - centroid-based.
- `--cluster-mode 1` connected component - transitive single linkage. **Use this.**
- `--cluster-mode 2` greedy incremental (CD-HIT-like) - centroid-based.

With centroid-based clustering, two sequences placed in *different* clusters can
still be similar to each other (each was merely closer to its own representative).
So whole-cluster assignment does not by itself keep splits apart. Connected-component
clustering uses transitive reachability (A~B, B~C => {A,B,C}), giving the strongest
guarantee. This is why the skill defaults to `--cluster-mode 1`.

## Tool choice

| Tool | Use it when | Notes |
|------|-------------|-------|
| **MMseqs2** `easy-cluster` | Default; below ~40% identity | "10000x faster than BLAST"; cascaded sensitive search captures distant homologs. |
| **MMseqs2** `easy-linclust` | Very large sets (10^6+) | Linear time; but practically bottoms out near ~50% identity. |
| **Graph-Part** | You want to *retain more data* than culling | Homology *partitioning*: keeps related sequences together while keeping as many as possible; guarantees separation at a chosen threshold. Good golden-reference comparator. |
| **CD-HIT** | Only above ~40% identity | Its k-mer short-word filter breaks below 40%; do not use it for sub-40% thresholds. |

Do **not** rely on thin pip MMseqs2 wrappers for reproducible splits - they are not
production-ready. Pin the official binary version (the skill records
`PINNED_MMSEQS_VERSION` in the provenance manifest); clustering results shift across
releases.

## Hermetic fallback

When no binary is on PATH, `scripts/cluster_sequence.py` falls back to k-mer Jaccard
similarity (`kmer_similar_pairs` -> `clusters_from_pairs`). This is monotonic with
identity and good enough to demonstrate and test the splitting logic, but it is an
approximation, not a substitute for MMseqs2 on real data. The report flags the
fallback loudly.

## Diagnostics to report

- Realized max train<->test identity and the full nearest-neighbour identity
  distribution (stratify metrics by it).
- A **BLAST/MMseqs2 nearest-neighbour transfer baseline** - the CAFA-style floor.
  A model that cannot beat it on a de-leaked split has shown nothing.
- Cluster-size distribution and the realized (often skewed) example-level ratio.

## Troubleshooting

| Problem | Fix |
|---------|-----|
| Everything collapses into one cluster | Threshold too high, or the data is genuinely a single family. Lower `--min-seq-id` or accept that no clean split exists; warn the user. |
| Realized test fraction much smaller than requested | A few giant clusters dominate. Expected with skewed family sizes; report the realized ratio. |
| Sub-40% threshold with CD-HIT gives odd clusters | CD-HIT cannot cluster below ~40%; switch to MMseqs2 cascaded or PSI-CD-HIT. |
| Passes identity split but downstream metrics still look inflated | Structural leakage below the identity threshold; add a Foldseek pass. |

## Citations

- Steinegger & Söding, MMseqs2, Nat. Biotechnol. 2017; Linclust, Nat. Commun. 2018.
- Rost, "Twilight zone of protein sequence alignments," Protein Eng. 12(2):85-94, 1999.
- Teufel et al., Graph-Part, NAR Genom. Bioinform. 2023, lqad088.
- Rao et al., TAPE, NeurIPS 2019 (25% identity splits).
