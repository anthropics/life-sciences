# Developability checks: what to compute and how (open tools only)

Everything here runs on CPU in seconds per design with python + biopython + freesasa (or biotite/prody) + DSSP, given the sequence and the designed (or re-predicted) complex model. Report values; compare to reference proteins or controls; do not invent pass/fail cut-offs. The sequence-only items below (motif counts and positions, pI and net charge, hydrophobic and charged windows, candidate protease sites, A280, cysteines) come out of `scripts/sequence_profile.py` as one table; exposure, interface membership and everything structure-based are computed as described here.

## 1. Chemical liabilities (sequence regex, then keep the exposed ones)
| Liability | Motif / rule | Notes |
|---|---|---|
| Asn deamidation | N-G (highest), N-S, N-T, N-N, N-H | rank NG ≫ NS > NT/NN; flexible, exposed loops deamidate fastest |
| Asp isomerisation | D-G (highest), D-S, D-T, D-D | iso-Asp formation; same exposure logic (the bundled script counts the two commonest, D-G and D-S) |
| Acid-labile bond | D-P | cleavage at low pH |
| Oxidation | M, W (also exposed H, C) with side-chain relative SASA > ~20 % | buried Met is low risk |
| N-glycosylation sequon | N-X-[S/T], X ≠ P | only matters if ever expressed in eukaryotic hosts; still report |
| Free cysteine | odd Cys count / unpaired Cys in the model | scrambling, dimerisation |
| N-terminal Gln / Glu | Q or E at position 1 | pyro-Glu formation |
| N-terminal Met processing | Met retained unless residue 2 is small (G, A, S, C, T, P, V) | affects the actual N-terminus of an E. coli product |
| Lys/Arg clusters, poly-ion runs, low complexity | ≥ 4 consecutive identical or same-charge residues; SEG-style low complexity | expression and non-specific binding flags |
For every hit record: position, relative SASA of the residue(s) (freesasa on the monomer model), secondary structure (DSSP), and whether it lies in the designed interface (any heavy atom within 4–5 Å of the target in the complex model). An NG in a buried helix core and an NG in an exposed interface loop are different risks — say which.

## 2. Charge
- pI and net charge at pH 7.4: biopython ProtParam (or Henderson–Hasselbalch with standard pKas). Report His count overall and at the interface (pH-sensitive binding; tumour and endosomal pH are acidic).
- Charge patches: assign +1 (K, R; H +0.1) / −1 (D, E) to exposed side chains (relative SASA > 25 %), cluster exposed charged atoms by distance (e.g. single-linkage at 6–8 Å), report the largest positive and largest negative patch as residue count and solvent-exposed area (Å²).

## 3. Hydrophobic surface, split by face
- Per-atom SASA with freesasa on the binder ALONE (from the complex model with the target deleted). Hydrophobic patch = connected set of exposed side-chain atoms of A, V, I, L, M, F, W, Y, P (C if present) within ~5 Å; report the largest patch area and total exposed hydrophobic area.
- SAP-like score per residue: sum over atoms within a 5 Å (or 10 Å) sphere of (atom SASA / max side-chain SASA of that residue type) × residue hydrophobicity (Black & Mould scale, Gly = 0 reference) — after Chennamsetty et al. 2009. Report the maximum per-residue SAP and list residues with SAP > 0 (a net-hydrophobic exposed neighbourhood, Chennamsetty's colouring convention) on the free surface; compare the count to the same number for your positive control / reference proteins rather than to a fixed cut-off.
- Do both calculations twice: for interface residues (within 4–5 Å of the target in the complex) and for the non-interface surface. The interface is allowed to be hydrophobic; a large hydrophobic patch on the free surface is the aggregation/non-specificity flag. Also report interface hydrophobic fraction (BindCraft keeps binder surface hydrophobicity < 35 %).

## 4. Aggregation-prone stretches (heuristic if no dedicated predictor is installed)
- Sliding-window mean hydrophobicity (Kyte–Doolittle, window 5–7) and β-propensity; flag windows that are both high and not buried in the model. If Aggrescan3D, CamSol, TANGO/WALTZ or NetSolP are installed (and licensed where required — TANGO, CamSol and NetMHCIIpan carry academic/commercial restrictions; Aggrescan3D standalone and NetSolP are openly distributed), run them and say which; otherwise say "heuristic only".

## 5. Protease exposure
Count canonical sites — trypsin (K/R not before P), chymotrypsin (F/W/Y/L), Asp-N (before D), Glu-C (after E) — that fall in DSSP loop/coil with relative SASA > 30 %. Hyperstable helical bundles resist proteolysis in practice; this is a flag list for the lab's protease/serum assay, not a prediction.

## 6. Stability and designability proxies
- Monomer re-prediction pLDDT and RMSD to the design (see thresholds.md).
- Inverse-folding log-likelihood / sequence recovery: score the designed sequence on the designed backbone with ProteinMPNN/SolubleMPNN (and an ESM-family LM pseudo-log-likelihood if available); re-design the backbone N times and report median sequence identity to the design and how many redesigns refold to the backbone — high recovery and many refolding alternatives indicate a designable, robust backbone.
- Sanity: helix capping residues present, no buried unsatisfied polar side chains in the core, no large internal cavities.

## 7. Format, termini, manufacturability
- Termini: relative SASA of the first/last 3 residues and their minimum distance to any target atom in the complex; recommend the terminus that is exposed and ≥ 15–20 Å from the interface for His/Strep tags, Fc, albumin-binding domains or PEG. If neither qualifies, say so — it changes how the protein can be made and half-life-extended.
- Host: Cys-free, no required glycosylation, < 15 kDa → E. coli cytoplasmic expression is the default assumption; note anything that argues otherwise.
- A280: count W and Y; if both are zero the protein cannot be quantified by A280 (common for small helical designs) — flag for the lab.
- Gene synthesis: for a codon-optimised coding sequence report GC % (aim 40–60 %), longest homopolymer (flag ≥ 8), direct/inverted repeats ≥ 12 bp — vendors' complexity screens reject on these.
- Intended format → PK note: state the intended use (bare binder, Fc fusion, albumin-binder fusion, multimer, imaging tracer) and the half-life class that implies (bare < 30–50 kDa: renal filtration, minutes–hours; albumin-binding: days; Fc: days–weeks). Do not claim a sequence-derived PK.

## 8. Immunogenicity (computable part only)
If an HLA class II binding predictor such as NetMHCIIpan, IEDB consensus or MixMHC2pred is installed (and licensed where required), scan 15-mers against a reference DRB1 panel and report predicted strong-binder density per 100 residues beside the same number for reference proteins run identically. Otherwise list under "not checked". Never describe a de novo sequence as "human-like"; the computable levers are epitope density and aggregation risk, and the decisive assays are in vitro (MAPPs, donor T-cell assays).
