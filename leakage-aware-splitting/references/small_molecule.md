# Small-molecule splitting (scaffold, Tanimoto, Butina)

## When to use

Molecular property / QSAR tasks: the model takes a molecule (SMILES) and predicts a
property. Random splits are far too optimistic here (Tox21 random ROC-AUC ~0.78 vs
scaffold ~0.68).

## The methods, from optimistic to realistic

1. **Bemis-Murcko scaffold split** (`scripts/split_small_molecule.py:scaffold_labels`).
   Reduce each molecule to its ring-and-linker core; assign whole scaffold groups to
   one split. The MoleculeNet/DeepChem baseline - better than random but still
   optimistic: benzene and pyridine are *different* scaffolds despite near-identical
   chemistry, and as scaffolds approach one-per-molecule, scaffold splitting converges
   to random. Use it as a baseline, never as the only control.
2. **Butina / fingerprint clustering** (`butina_labels`) on ECFP4 Morgan fingerprints,
   then whole-cluster assignment. Harder and more realistic; the cutoff is a distance
   (1 - Tanimoto). O(n^2), fine at benchmark scale.
3. **Tanimoto nearest-neighbour control** (`tanimoto`). Always report, for each test
   molecule, the max ECFP4 Tanimoto to any train molecule. The Lo-Hi splitter enforces
   no test molecule above 0.4 Tanimoto to train (motivated by an EMA expert study where
   ~half of medicinal chemists call a pair "different" below 0.4 Tanimoto).

## Edge cases the skill handles

- **Acyclic molecules** have an empty Bemis-Murcko scaffold (RDKit returns ""). They
  get their own singleton clusters, never pooled together - pooling them would
  manufacture false similarity.
- **Invalid SMILES** also return "" and are isolated as singletons rather than crashing.

## Decoy / analog bias (virtual screening)

If the dataset has actives + decoys, scaffold disjointness is not enough. DUD-E rewards
analog bias (a few chemotypes dominate actives) and decoy bias (decoys trivially
separable in 2D) - CNNs hit AUC > 0.9 by exploiting these, not by learning physics.
Report the **AVE bias** metric (`leakage_metrics.ave_bias`) and prefer AVE-debiased /
LIT-PCBA benchmarks over DUD-E. Caveat (Sundar & Colwell): over-aggressive AVE
debiasing that deletes near-neighbours can *worsen* far-test generalization, so treat a
high-AVE warning as "investigate", not "always split harder".

## Citations

- Bemis & Murcko, J. Med. Chem. 1996; Wu et al., MoleculeNet, Chem. Sci. 2018.
- Steshin, Lo-Hi splitter, NeurIPS 2023 (0.4 Tanimoto).
- Wallach & Heifets, AVE bias, JCIM 2018; Chen et al., DUD-E hidden bias, PLoS One 2019.
- Sundar & Colwell, debiasing can hurt generalization, JCIM 2020.
