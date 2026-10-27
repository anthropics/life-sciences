# Dossier template (a markdown dossier plus a per-design CSV)

Write for a reviewer who did not watch you work. Numbers for controls sit beside numbers for designs. Every section ends with what it does NOT establish.

## 0. Summary (5–8 lines)
How many designs pass your pre-stated consensus rule; the controls' scores; how many candidates you set aside and the commonest reason; the two or three developability flags that recur; what was not checked.

## 1. Binding evidence beyond the objective
Predictors used in design/scoring vs predictors used here (table). Per design: interface confidence per predictor × seed summary, DockQ / iRMSD to the designed pose, consensus call. Controls through the identical protocol: positive(s) (name, provenance, e.g. "chain B of PDB xxxx"), scrambles, unrelated folded protein, decoy target — with their numbers. One sentence on independence (which of these models share training data or architecture with your design tools).

## 2. Function
Superposition test against the natural partner/known ligand complex: clash volume or clashing-atom count, footprint overlap fraction, "predicted blocker yes/no" per design.

## 3. Specificity and species cross-reactivity
Per design × per panel protein: confidence and whether a pose forms at all. Paralog/off-target rows: want low. Ortholog rows: report as predicted cross-reactivity (wanted or not, per the use case). State that co-folding confidence separates on- from off-target only roughly.

## 4. Robustness
Per design: pose persistence across seeds/models (contact-persistence fraction, pose RMSD spread); self-association co-fold result; monomer refold pLDDT/RMSD. For the leads: proposed binding-dead mutant(s) and their re-predicted scores.

## 5. Native context
Target glycosylation sequons/annotated sites, disulfides, natural variants within or adjacent to each footprint; result of re-predicting leads against the full extracellular region / full domain architecture; approach direction relative to membrane or neighbouring domains.

## 6. Developability (per design table + prose on recurring flags)
Liabilities (exposed, with interface membership), pI, net charge, largest ± charge patches, hydrophobic patch and SAP-like maxima for interface vs free surface, aggregation-prone windows, protease-exposed sites, stability/designability proxies, free terminus recommendation, host/A280/gene-synthesis notes, interface His, intended format and PK implication, epitope density if computed. Compare to reference proteins if you have them; no pass/fail language without a sourced threshold.

## 7. Interface metrics (uncalibrated unless Rosetta/FoldX available)
Buried SASA, shape complementarity, H-bonds/salt bridges, buried unsatisfied polar atoms, interface hydrophobic fraction — designs and positive control side by side.

## 8. Candidates set aside
Table: candidate, the check it failed, the number, what replaced it.

## 9. Not checked, and why
Every item from Parts A and B of the skill that you did not run (tool absent, licence, time, no sequence available), and every lab-only property (expression/SEC, affinity/kinetics, specificity panel, Tm, protease/serum stability, structure, immunogenicity assays, PK).

## The per-design CSV
One row per design; first column the binder sequence (so rows can be joined to whatever table lists the designs); then one column per reported quantity with self-explanatory names and units in the header (e.g. af2ig_ipae_mean, boltz_dockq_to_design_min, monomer_plddt, monomer_rmsd_A, consensus_pass, blocker, homodimer_iptm, n_exposed_NG, pI, net_charge_pH7.4, max_pos_patch_A2, max_hphob_patch_free_A2, sap_max_free, free_terminus, a280_measurable, …). Include the controls as extra rows with a 'role' column (design / positive_control / scramble / decoy).
