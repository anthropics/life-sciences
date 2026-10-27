# Back-tests — what the screen sees, and a public reference panel for developability

> Everything on this page is computed: structure-predictor confidence, and descriptors calculated from sequences and structures, for public reference proteins only. Nothing here is a binding, expression, stability or safety measurement, no reference protein was assayed for this page, and no design from any session appears on it: results from real runs live in `example-run/` and `evals/`, never here.

Two questions arrive with every panel. Does the screen that scored it recognise proteins known to bind this target, and reject proteins that do not? And on the properties a developability group computes first, do the designs resemble proteins that have been given to patients? The rubric's Part 3 answers the first with four controls; section 1 widens it to marketed antibodies, published de novo binders and format-matched decoys, all public sequences. The second has no answer key for miniproteins, so section 2 gives a public reference panel of marketed, clinical-stage and research small proteins with their computed descriptor ranges, for your dossier to be read against. Both were computed on the public PD-L1 sample target with open tools and no design session in the loop. Read them as a description of the instrument and a yardstick, not as evidence about any panel.

## 1. Known binders and non-binders through the same screen

### Method

Twenty-six public sequences went through the screen exactly as `rubric.md` Part 3 describes it: ESMFold2 in its fast single-sequence configuration, ESMFold2 in its full configuration with the target's alignment, and Protenix with the target's alignment, against the same 115-residue PD-L1 construct (UniProt Q9NZQ7, residues 18–132); five seeds per program, best seed kept; interface confidence is ipSAE (Dunbrack, 2025), the lower of the two chain directions; "all three" is the mean of the three programs. The three programs were run in accelerated open builds of the kind the guide's section 4 links; other builds of the same programs may give values that differ in the second decimal. The binder side is always folded from its sequence alone, and no similarity filter was applied. Every sequence comes from a public entry:

- **The Part 3 anchors**, measured again in the same batch: natural PD-1 (PDB 4ZQK, chain B), the engineered high-affinity PD-1 (PDB 5IUS, chain A), shuffled natural PD-1, barstar (PDB 1BRS, chain D), and barstar with its own partner barnase.
- **Known PD-L1 binders in other formats**: the single-domain antibody KN035 (PDB 5JDS, chain B, tags removed), which is the PD-L1-binding domain of envafolimab; the variable regions (Fv) of atezolizumab, durvalumab, avelumab and BMS-936559, taken from their PD-L1 co-crystal entries (PDB 5X8L, 5X8M, 5GRJ, 5GGT) and checked against the INN sequence where ChEMBL holds one (CHEMBL3707227, CHEMBL3301587); and the published de novo PD-L1 binders DBL1_03 (PDB 7XYQ, chain B) and DBL2_02 (PDB 7XAD, chain C) of Gainza et al., *Nature* 2023.
- **Presumed non-binders** matched to those formats. Binders of unrelated targets: one single domain of caplacizumab (PDB 7EOW, chain B), the single-chain Fv brolucizumab (INN, CHEMBL3707357), an anti-EGFR adnectin (PDB 3QWQ, chain B), the anti-HER2 DARPin H10-2-G (PDB 4HRN, chain A), the anti-HER2 affibody ZHER2:342 (PDB 3MZW, chain B), an anti-CTLA-4 anticalin (PDB 3BX7, chain A) and the de novo interleukin-2 mimetic Neoleukin-2/15 (PDB 6DG6). Folded proteins with no binding role: Top7 (PDB 1QYS), α3D (PDB 2A3D), protein G B1 domain (PDB 1PGA), ubiquitin (PDB 1UBQ). None has, to our knowledge, been assayed against PD-L1; "non-binder" is a presumption.

The screen takes one chain per binder, so each Fv was made single-chain twice: heavy variable domain, a (Gly₄Ser)₃ linker, light variable domain; and the reverse. These constructs are the guide's, made for this page. None is a marketed format, no linker is resolved in any of the deposited structures, and the two orientations of one antibody are not independent observations. Tags were removed where the entry marks them; DBL1_03 keeps a five-residue N-terminal remnant, per the depositors' numbering.

### Results

| Role | Protein (public source) | Residues | ESMFold2, fast | ESMFold2, full + alignment | Protenix + alignment | All three |
|---|---|---|---|---|---|---|
| anchor: binder | Natural PD-1 ectodomain (PDB 4ZQK, chain B) | 118 | 0.81 | 0.82 | 0.85 \* | 0.83 |
| anchor: binder | Engineered high-affinity PD-1 (PDB 5IUS, chain A) | 129 | 0.82 | 0.82 | 0.87 | 0.84 |
| anchor: negative | Shuffled natural PD-1 | 118 | 0.00 | 0.00 | 0.00 | 0.00 |
| anchor: negative | Barstar (PDB 1BRS, chain D), against PD-L1 | 89 | 0.00 | 0.01 | 0.07 | 0.03 |
| known binder | KN035 single-domain antibody (PDB 5JDS, chain B) | 128 | 0.02 | 0.28 | 0.01 | 0.10 |
| known binder | Atezolizumab Fv, heavy–linker–light (PDB 5X8L; INN) | 240 | 0.10 | 0.02 | 0.02 | 0.05 |
| known binder | Atezolizumab Fv, light–linker–heavy | 240 | 0.16 | 0.01 | 0.02 | 0.07 |
| known binder | Durvalumab Fv, heavy–linker–light (PDB 5X8M; INN) | 244 | 0.01 | 0.00 | 0.22 | 0.08 |
| known binder | Durvalumab Fv, light–linker–heavy | 244 | 0.01 | 0.00 | 0.01 | 0.01 |
| known binder | Avelumab Fv, heavy–linker–light (PDB 5GRJ) | 245 | 0.81 | 0.78 | 0.88 \* | 0.83 |
| known binder | Avelumab Fv, light–linker–heavy | 245 | 0.79 | 0.79 | 0.88 \* | 0.82 |
| known binder | BMS-936559 Fv, heavy–linker–light (PDB 5GGT) | 244 | 0.78 | 0.68 | 0.46 | 0.64 |
| known binder | BMS-936559 Fv, light–linker–heavy | 244 | 0.78 | 0.68 | 0.58 \* | 0.68 |
| known binder | DBL1_03, de novo (PDB 7XYQ, chain B) | 124 | 0.02 | 0.66 | 0.33 | 0.34 |
| known binder | DBL2_02, de novo (PDB 7XAD, chain C) | 105 | 0.81 | 0.78 | 0.84 | 0.81 |
| presumed non-binder | Top7, de novo fold (PDB 1QYS) | 100 | 0.17 | 0.11 | 0.44 | 0.24 |
| presumed non-binder | α3D, de novo three-helix bundle (PDB 2A3D) | 73 | 0.03 | 0.00 | 0.01 | 0.01 |
| presumed non-binder | Protein G, B1 domain (PDB 1PGA) | 56 | 0.00 | 0.15 | 0.02 | 0.05 |
| presumed non-binder | Ubiquitin (PDB 1UBQ) | 76 | 0.01 | 0.01 | 0.13 | 0.05 |
| presumed non-binder | Neoleukin-2/15, de novo interleukin-2 mimetic (PDB 6DG6, chain A) | 100 | 0.00 | 0.01 | 0.02 | 0.01 |
| presumed non-binder | Caplacizumab, one single domain (PDB 7EOW, chain B) | 128 | 0.00 | 0.00 | 0.03 | 0.01 |
| presumed non-binder | Brolucizumab, single-chain Fv (INN, CHEMBL3707357) | 252 | 0.01 | 0.01 | 0.42 | 0.15 |
| presumed non-binder | Anti-EGFR adnectin (PDB 3QWQ, chain B) | 108 | 0.22 | 0.12 | 0.33 | 0.22 |
| presumed non-binder | Anti-HER2 DARPin H10-2-G (PDB 4HRN, chain A) | 124 | 0.01 | 0.22 | 0.25 | 0.16 |
| presumed non-binder | Anti-HER2 affibody ZHER2:342 (PDB 3MZW, chain B) | 58 | 0.02 | 0.05 | 0.13 | 0.06 |
| presumed non-binder | Anti-CTLA-4 anticalin (PDB 3BX7, chain A) | 178 | 0.00 | 0.00 | 0.01 | 0.00 |
| scale | Barstar with its own partner barnase (PDB 1BRS), a different target | 89 + 110 | 0.89 | 0.89 † | not run | — |

\* Protenix's best-of-five value rests on one or two seeds; the other seeds gave near-zero interface confidence (natural PD-1: one seed at 0.85, four at 0.01; avelumab: two seeds and one seed in the two orientations; BMS-936559 light–linker–heavy: one seed). Both ESMFold2 configurations placed every binder they recognised on all five seeds. † Single-sequence mode: no alignment exists for barnase, and Protenix was not run on this pair for the same reason. The anchors reproduce the Part 3 table, measured the day before, to within 0.02 in every entry, so the two batches are on one ruler.

### Reading

1. **Non-binders stay low.** The screen kept all thirteen negatives (the eleven presumed non-binders and the two negative anchors) below 0.25 on "all three" (highest: Top7 0.24, the adnectin 0.22, the DARPin 0.16, brolucizumab 0.15), and none reaches 0.5 on any seed of any program. The two antibody-format negatives (0.01, 0.15) score like the rest, so antibody format alone does not produce false positives. Protenix on its own is the noisiest on negatives (best seeds 0.44, 0.42, 0.33): the clean separation belongs to the three-program mean, not to any single program.
2. **Known binders split by format.** Counting each antibody once, five of nine score with the binder anchors: natural and engineered PD-1 (0.83, 0.84), avelumab's Fv (0.83, 0.82), DBL2_02 (0.81) and BMS-936559's Fv (0.64, 0.68). One is marginal: DBL1_03 at 0.34, placed by the full ESMFold2 configuration on all five seeds (0.57–0.66) and by neither other program. Three score like non-binders: KN035 (0.10), atezolizumab's Fv (0.05, 0.07) and durvalumab's Fv (0.08, 0.01). Six of nine therefore sit above every negative, DBL1_03 by only 0.10; the misses are the Fvs of two marketed anti-PD-L1 antibodies and the binding domain of a third, envafolimab, approved in China.
3. **The misses look like docking misses, not folding misses.** The missed single-chain Fvs are predicted with about the chain confidence of the recognised ones (binder pLDDT 0.86–0.89 against 0.89–0.91 in the fast configuration — a chain average dominated by framework β-sheet, silent about the loops); KN035 lower (0.77–0.79; 0.50–0.64 on Protenix). They are simply not placed on PD-L1: best-seed ipTM 0.19–0.69 across programs, ipSAE never above 0.28. This is a blind spot of the screen for antibody-loop interfaces — patchy rather than uniform, since avelumab and BMS-936559, built the same way from entries of the same period, are found — and not news about the antibodies. Nothing at sequence level here explains which Fvs are found.
4. **Best-of-five hides how marginal some hits are.** Natural PD-1's Protenix value rests on one seed of five; avelumab's on one or two; BMS-936559's heavy–linker–light construct on none (best 0.46). DBL2_02 and the engineered PD-1 are placed on all five seeds by all three programs.
5. **Confidence, not affinity.** Where a binder is recognised the value runs 0.64–0.84 and does not track binding strength: the engineered PD-1, reported to bind thousands of times tighter than natural PD-1, scores 0.84 against 0.83; the one recognised marketed antibody 0.83; the antibody that did not go beyond early-phase trials, 0.64–0.68. Expect a panel of designed helical miniproteins to score *above* every real binder in this table, and do not read that as tighter binding: idealised de novo proteins folded from one sequence are what these programs are most confident about, natural domains and antibody loops are what they are least sure of, and a panel selected with these same programs scores well on them by construction.
6. **So use it as a floor-and-ceiling check.** A design that scores with the negatives has failed. A design that scores with the binder anchors has earned a look by a protein scientist, nothing more. Do not rank designs against one another by "strength" on it, and never against antibodies: it missed three of the five antibody-derived binders tested here, three of the four marketed ones.

### Caveats that travel with these numbers

- Small numbers: nine binders in four formats, thirteen negatives, no error bars; the two orientations of each Fv are one observation.
- The single-chain Fv constructs are the guide's, not a marketed format; KN035 is the camelid-framework domain in PDB 5JDS, and the domain in marketed envafolimab may differ.
- Every binder here has a public PD-L1 co-crystal structure, and structure predictors are trained on that public database, so some may be partly recalled. Recall did not rescue atezolizumab, durvalumab or KN035; it may flatter the binders that were found.
- "Recognised" means a confident interface was output. We did not compare the predicted avelumab, BMS-936559 or DBL2_02 complexes with PDB 5GRJ, 5GGT or 7XAD; a DockQ calculation would say whether the confident poses are the right ones.
- Folding the binder without an alignment favours idealised de novo proteins over natural domains and antibodies — here and in every panel scored this way.
- The negatives are presumed, not assayed.

## 2. Developability descriptors: a public reference panel

### Why a panel, not a threshold

There is no validated computational developability rule set for de novo miniproteins. Antibody rule sets (the Therapeutic Antibody Profiler, clinical-stage assay panels) are calibrated on IgGs and do not transfer to 6–13 kDa helical bundles. What a reviewer can do is compute the descriptors a developability group looks at first, for the designs *and* for small proteins that have been given to patients, and read one against the other. This section supplies the second half: a public reference panel and its computed ranges. The `binder-gauntlet` skill's Part B computes the same descriptors for your designs and reports them beside these (its Principle 4); nothing here is a pass line.

### The reference proteins

Twenty-four proteins with public sequences: sixteen marketed or clinical-stage, and eight research proteins (published de novo binders and folds, and engineered scaffolds standing in for clinical scaffold drugs of the same classes whose sequences are not public). "Status" is public record from labels and company statements, given so the table can be read; it was not re-verified for this page. Failures rarely publish sequences, and the failures and safety cases that could be included failed for reasons no sequence descriptor measures (pre-existing antibodies to a single-chain Fv framework and an intravitreal route; efficacy; portfolio decisions).

| Reference protein | Status | Sequence | Structure used (partner) |
|---|---|---|---|
| *Marketed or clinical (16)* | | | |
| Anakinra (IL-1 receptor antagonist) | marketed | INN (CHEMBL1201570) | PDB 1IRA (IL-1 receptor type 1) |
| Aprotinin | marketed; suspended for some years after a safety trial | PDB 4PTI = INN (CHEMBL1201619) | PDB 4PTI (none) |
| Desirudin (hirudin variant 1) | marketed; later commercially discontinued | PDB 2HIR = ChEMBL hirudin record | PDB 4HTC (thrombin); a hirudin variant differing at 5 of 61 resolved positions |
| Ecallantide (engineered Kunitz domain) | marketed; anaphylaxis warning on its label | INN (CHEMBL1201837) | none |
| Caplacizumab, bivalent | marketed | INN (CHEMBL2109624) | none |
| Caplacizumab, one single domain | domain of the above | PDB 7EOW = INN unit | PDB 7EOW (von Willebrand factor A1) |
| Ozoralizumab, trivalent single-domain | marketed (Japan) | INN (CHEMBL1743054) | none |
| Vobarilizumab, bispecific single-domain | stopped after phase 2 | INN (CHEMBL3833343) | none |
| Brolucizumab, single-chain Fv | marketed; post-marketing immunogenicity and retinal-vasculitis signal | INN (CHEMBL3707357) | none |
| KN035 single-domain antibody | binding domain of envafolimab (approved in China) | PDB 5JDS, tags removed | PDB 5JDS (PD-L1) |
| Atezolizumab Fv (single-chain construct, as in section 1) | Fv of a marketed antibody | INN (CHEMBL3707227) = PDB 5X8L | PDB 5X8L (PD-L1) |
| Durvalumab Fv (construct) | Fv of a marketed antibody | INN (CHEMBL3301587) = PDB 5X8M | PDB 5X8M (PD-L1) |
| Avelumab Fv (construct) | Fv of a marketed antibody | PDB 5GRJ | PDB 5GRJ (PD-L1) |
| BMS-936559 Fv (construct) | Fv of an antibody that did not go beyond early-phase trials | PDB 5GGT | PDB 5GGT (PD-L1) |
| Neoleukin-2/15 (de novo interleukin-2 mimetic) | parent of a clinical candidate discontinued in phase 1 | PDB 6DG6, chain A | PDB 6DG5 (IL-2 receptor β and γ chains) |
| Taldefgrobep alfa, adnectin domain only | domain of an Fc fusion whose phase 3 missed its primary endpoint | INN (CHEMBL4297914), residues 244–340 | none |
| *Research (8)* | | | |
| DBL1_03, DBL2_02 (de novo PD-L1 binders) | research | PDB 7XYQ, 7XAD | same entries (PD-L1) |
| Top7 (de novo fold, not a binder) | research | PDB 1QYS | PDB 1QYS (none) |
| Anti-EGFR adnectin | research; class of a former clinical adnectin | PDB 3QWQ | PDB 3QWQ (EGFR) |
| Anti-HER2 DARPin H10-2-G | research; class of clinical DARPins | PDB 4HRN | PDB 4HRN (HER2 domain IV) |
| Consensus DARPin (unselected) | research | PDB 2QYJ | PDB 2QYJ (none) |
| Anti-HER2 affibody ZHER2:342 | research; lineage of a clinical imaging agent | PDB 3MZW | PDB 3MZW (HER2) |
| Anti-CTLA-4 anticalin | research | PDB 3BX7 | PDB 3BX7 (CTLA-4) |

Not size-matched and so left out of every range, though a reviewer may want them as context: the full heavy and light chains of atezolizumab and durvalumab, the full 340-residue taldefgrobep alfa chain (half of it an antibody constant region), and obligate-oligomer subunits such as designed nanoparticle components, whose free surface is an assembly interface.

### Descriptors

Computed with Biopython from the sequence and, where one exists, from a structure (Shrake–Rupley accessibility, 1.4 Å probe, protein heavy atoms only). "Exposed" means relative accessibility ≥ 0.20 in the binder alone; an "interface residue" loses at least 1 Å² when the partner is present; a "surface region" is, for each exposed residue of a class, the summed accessible surface of all exposed residues of that class whose Cβ lies within 10 Å of it, maximised over residues — a residue-level analogue of the spatial-aggregation-propensity neighbourhood sum, not a clustered patch. The skill's `scripts/sequence_profile.py` reproduces every sequence column; the structure columns need a solvent-accessibility library and a model of the complex.

| Group | Descriptors |
|---|---|
| Charge | isoelectric point (pI); net charge at pH 7.4; fraction of charged residues (D, E, K, R); largest positive (K, R) and negative (D, E) surface regions; extreme net charge over any 10-residue window; longest charged run |
| Sequence hydrophobicity | GRAVY (mean Kyte–Doolittle hydropathy); longest run of A/F/I/L/M/V/W/Y; highest mean hydropathy over any 7-residue window; number of 7-residue windows with mean ≥ 1.5; aromatic fraction (F, W, Y) |
| Hydrophobic surface | share of the unbound surface carried by F/I/L/M/V/W/Y residues (alanine deliberately excluded); largest hydrophobic surface region, unbound, and the same with interface residues removed; hydrophobic share of the surface left exposed in the complex |
| Chemical liabilities | counts of Asn-Gly, Asn-Ser, Asp-Gly, Asp-Ser, Asp-Pro, N-X-S/T sequons and N-terminal Gln, and their total; Asn-Gly/Asn-Ser/Asp-Gly/Asp-Ser sites that are exposed, and that sit in the interface; Met and Trp exposed, and in the interface; cysteine count and cysteines unpaired in the structure |
| Size and interface | length; molecular weight; surface buried in the interface, binder side |
| Other | lysine count; absorbance at 280 nm of a 1 mg/ml solution from Trp and Tyr content (zero means not measurable); helix-favouring residue share; instability index (shown because it is asked for, not because it discriminates) |

**Structures are not like for like with a design's.** Reference structure columns come from experimental structures, which lack disordered loops and termini (a median 6% of each sequence unresolved) and so read slightly more hydrophobic than complete models; a design's come from a predicted complex. Eighteen of the 24 references have a structure and 15 a partner chain; antibody Fab chains were cut to the Fv. Say which is which whenever the two sit in one table.

### Reference ranges

| Descriptor | All 24 references: n | range | The 16 marketed or clinical: n | range |
|---|---|---|---|---|
| Length (residues) | 24 | 58 to 363 | 16 | 58 to 363 |
| Molecular weight (kDa) | 24 | 6.5 to 38.4 | 16 | 6.5 to 38.4 |
| Isoelectric point | 24 | 4.04 to 9.24 | 16 | 4.04 to 9.24 |
| Net charge at pH 7.4 | 24 | −12.2 to 5.54 | 16 | −9.59 to 5.54 |
| Fraction of charged residues | 24 | 0.13 to 0.34 | 16 | 0.13 to 0.34 |
| Longest charged run | 24 | 2 to 5 | 16 | 2 to 5 |
| Largest positive region (Å²) | 18 | 145 to 687 | 10 | 214 to 445 |
| Largest negative region (Å²) | 18 | 172 to 582 | 10 | 173 to 582 |
| GRAVY | 24 | −0.90 to 0.04 | 16 | −0.77 to −0.11 |
| Longest hydrophobic run | 24 | 3 to 6 | 16 | 3 to 6 |
| Highest 7-residue hydropathy | 24 | 1.11 to 2.74 | 16 | 1.11 to 2.74 |
| Hydrophobic 7-residue windows (count) | 24 | 0 to 12 | 16 | 0 to 11 |
| Aromatic fraction (F, W, Y) | 24 | 0.032 to 0.150 | 16 | 0.046 to 0.150 |
| Absorbance at 280 nm, 1 mg/ml | 24 | 0.23 to 2.33 | 16 | 0.43 to 2.33 |
| Hydrophobic share of unbound surface | 18 | 0.151 to 0.271 | 10 | 0.151 to 0.234 |
| Largest hydrophobic region outside the interface (Å²) | 15 | 53 to 353 | 9 | 102 to 283 |
| Hydrophobic share of bound surface | 15 | 0.10 to 0.23 | 9 | 0.12 to 0.18 |
| Interface buried on the binder side (Å²) | 15 | 572 to 1,744 | 9 | 590 to 1,744 |
| Liability motifs, total | 24 | 0 to 9 | 16 | 0 to 9 |
| Helix-favouring residue share | 24 | 0.167 to 0.500 | 16 | 0.192 to 0.500 |
| Instability index | 24 | 7.5 to 65.1 | 16 | 11.5 to 65.1 |
| Lysine count | 24 | 2 to 14 | 16 | 2 to 11 |

### Reading a dossier against it

- **"Outside" is a prompt to look, not a verdict**, and on several rows outside is the favourable side (less hydrophobic surface, fewer liabilities). Sixteen survivors validate no threshold.
- **Charge.** Nine of the ten antibody-derived references are basic (pI 7.9–9.2; brolucizumab, 4.9, is the exception), while anakinra, ecallantide, desirudin and every de novo or scaffold reference are acidic (4.04–6.75). A helical de novo panel will usually sit with the acidic group; what deserves a look is a design beyond every reference on fraction of charged residues or on the size of a charged surface region, since large charged patches are a recognised non-specificity and clearance flag for antibodies and are untested for this format.
- **Absorbance at 280 nm.** Every reference has at least one Trp or Tyr (the weakest, the anti-HER2 DARPin, absorbs 0.23 at 1 mg/ml). De novo helical designs often have neither, and then routine A280 concentration measurement does not work. It is trivial to design out at round one and a nuisance to discover at expression.
- **Labile Asn/Asp motifs in the interface.** Across the 15 references with a complex structure, 4 of their 31 Asn-Gly, Asn-Ser, Asp-Gly and Asp-Ser sites (13%) sit in the binding interface. A design that puts such a site in its predicted interface has a modification risk that touches binding; those are the sites to stress-test first, or to design out while it is cheap.
- **Sequence hydrophobicity and alanine.** The Kyte–Doolittle scale counts alanine as hydrophobic; the surface measures above do not. Helical de novo designs are often alanine-rich, so a design can exceed every reference on GRAVY or on the 7-residue window while sitting inside the range on hydrophobic surface. Report both, say which residues drive the sequence number, and let a solubility assay or a dedicated aggregation predictor settle it.
- **Instability index** exceeds its textbook cut-off of 40 for ten of the sixteen marketed or clinical references. Show it because reviewers ask; do not act on it.

> **What this comparison cannot say.** It places a design inside, or flags where it sits outside, the range of computed descriptors spanned by small proteins that reached patients. Sixteen survivors validate no threshold, and there is no validated computational developability rule set for de novo miniproteins. De novo designs look clean on liabilities and charge balance partly by construction; that is real, and it is not "more developable than anakinra". Solubility, expression, aggregation, stability and protease resistance are measured properties: descriptors flag risks and cannot clear them. Pharmacokinetics at this size is set by format, not sequence — an unmodified protein this small is filtered by the kidney within hours, and half-life comes from an Fc, an albumin binder or a polymer — so none is predicted. Immunogenicity has a computable part (T-cell epitope content, aggregation) that is not computed here, and decisive parts (impurities, route, target biology, pre-existing antibodies) that are not computable. A design's surface numbers come from a predicted complex; the references' from experimental ones.

**Not computed here, and why.** Dedicated aggregation and solubility predictors — Aggrescan3D, NetSolP, CamSol and TANGO — carry terms a customer must check before sending proprietary sequences (CamSol is a web server whose terms exclude them; TANGO carries an academic licence); the sequence-window heuristic stands in until they are in your `TOOLS.md`. T-cell epitope load (NetMHCIIpan, academic licence). Protease-site exposure and anything needing secondary-structure assignment; designability and inverse-folding likelihoods. A calibration of any descriptor against measured expression or solubility: public datasets of laboratory-tested de novo designs exist (Proteinbase is one) and a customer with the right to use one could compute it (check the dataset's terms for your use first); none is reproduced here. Anything pharmacokinetic, for the reason in the box.

## 3. How this maps onto the binder-gauntlet skill

The guide ships the `binder-gauntlet` skill and tells Claude to run its cheap checks at every selection round, not only at the end. Section 1 is the skill's **Part A3, controls through the identical protocol**, at larger scale than a session would afford (one A3 item, a decoy *target* of PD-L1's fold class, is not covered here); section 2 is the reference half of its **Part B, developability profile** — B1 liabilities with exposure and interface membership, B2 charge, B3 hydrophobic surface by face, the heuristic half of B4 and the A280 and host notes of B8 — against which the skill reports a design's own numbers, as its Principle 4 asks (B5–B7 and B9 need tools this page does not assume). The three items in "Reading a dossier against it" that recur in helical de novo panels — no Trp or Tyr, labile Asn/Asp motifs in the interface, alanine-heavy stretches that read hydrophobic on sequence scales — cost a redesign when found at the end and nothing when found at round one, which is why the guide puts them inside the loop.

## 4. Reproducing this on your own target

- For section 1, swap in your target's sequences: its natural partner from a deposited complex, any public engineered binder, known binders in whatever formats exist (antibodies made single-chain the same way, labelled as your construction), a shuffle, two or three unrelated folded proteins of similar size, and same-format binders of *other* targets as decoys.
- Run them through the screen your session used for its panel, unchanged — programs, settings, construct, alignment, seeds, best-of-n rule — and report per program and combined beside your panel's range, in the Part 3 table format. The 26 public sequences of section 1 took under an hour on one data-centre GPU.
- For section 2, keep the panel above if your designs are small proteins, or pick reference proteins of your own modality and size with public sequences (INN records in ChEMBL; deposited structures with the partner where they exist); run the skill's Part B on them and on your designs, and report ranges side by side, saying which structure columns are predicted and which experimental.
- Keep the limits attached when the numbers travel: presumed negatives, constructs that are not marketed formats, survivors that validate no threshold, confidence that is not affinity.
