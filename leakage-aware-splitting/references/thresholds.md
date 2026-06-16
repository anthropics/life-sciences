# Similarity thresholds: conventions, not laws

Every threshold below is a widely used heuristic grounded in a specific finding.
None is a law. The "right" cutoff depends on what generalization you are testing;
expose it as a parameter and document the choice in the provenance manifest.

| Modality | Threshold | What it means | Grounding |
|----------|-----------|---------------|-----------|
| Protein sequence | **30% identity** (often 25%) | Below this, sequence comparison can no longer reliably establish homology | Rost 1999 twilight zone: ">~30% identity, 90% of pairs homologous; <25%, <10% are." TAPE uses 25%. Structure-grade work pushes to <20%. |
| Nucleotide sequence | **higher than protein** | Random identity baseline is ~25% for a 4-letter alphabet (vs ~5% for 20 amino acids) | Set well above the random floor; Graph-Part supports DNA/RNA. |
| Protein structure | **TM-score >= 0.5** | Above 0.5 two structures share the same fold | Standard fold-equivalence threshold; Foldseek also uses E-value < 0.01. |
| Small molecule | **ECFP4 Tanimoto < 0.4** (test vs train) | Below ~0.4 medicinal chemists consider two molecules "different enough" | Lo-Hi splitter (Steshin, NeurIPS 2023), replicating an EMA expert study. |
| Small molecule (clustering) | **Tanimoto ~0.6** Butina cutoff | Cluster cutoff that empirically works for Butina splits | DeepChem ButinaSplitter. |
| Redundancy reduction (mild) | 50% / 70% / 90% | Removing near-duplicates, NOT testing generalization | PDB standard clusterings at 30/50/70/90/95/100%. |

## Why "necessary but not sufficient"

A sequence-identity split can still leak structurally: two proteins below 30%
identity can share a fold or binding interface (the PPI-leakage example: PDB 1K3F
and 1K9S share only 26.5% identity yet are structural near-duplicates). Scaffold
disjointness can still leave high-Tanimoto pairs across splits (benzene vs pyridine
are "different scaffolds" but near-identical chemistry). **Combine metrics per
modality**; never rely on a single axis.

## Picking a threshold

- Match the split to the **intended deployment**. If the model will be used on
  same-family inputs, an aggressive low-identity split understates true utility.
- Stricter is not automatically more honest (Sundar & Colwell): over-aggressive
  debiasing that deletes near-neighbours can *reduce* far-test performance. Validate
  that a split predicts prospective performance, not merely that it lowers scores.
- When unsure, report metrics **stratified by nearest-train similarity bins** so the
  reader can see performance as a function of distance rather than trusting one
  headline number.
