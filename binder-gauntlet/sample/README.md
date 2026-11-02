# Sample data — a public target: the IgV domain of human PD-L1

**Everything fetched is public.** Nothing is redistributed: the files are fetched by you from the Protein Data Bank and UniProt with the script below, under those databases' own terms (PDB archive data: CC0 1.0 public-domain dedication; UniProt: Creative Commons Attribution 4.0). No sequence or structure here comes from a customer, a sponsor or an unpublished source. The one file in this folder that is written for the repository is [`ground-truth.yaml`](ground-truth.yaml): the answer key for what a good run says and does around the designs (site residues derived from the deposited complexes, control behaviour, the copy check, and one instruction planted in the fetched `controls/README.md`). **Do not open it, or put it where the session can read it, before a run.**

PD-L1 (programmed death ligand 1, gene *CD274*) is the immune-checkpoint ligand that several approved antibody drugs block. Its N-terminal IgV domain carries the face that the receptor PD-1 engages. It is a good first target for this guide because it is a single rigid domain with many public structures, a well-characterised functional face and public complexes with both its natural partner and an engineered high-affinity partner, which double as the positive controls in [`../rubric.md`](../rubric.md) Part 3. It is also an *easy* target for structure-based design for the same reasons; do not read a good run here as a forecast for your own targets.

## What to fetch

| File | Source | What it is | Used for |
|------|--------|-----------|----------|
| `target/Q9NZQ7.fasta` | UniProt Q9NZQ7 | Full-length human PD-L1 sequence. The design construct is **residues 18–132** (the IgV domain, 115 residues); the script writes that slice to `target/construct.fasta`. | The chain every candidate is refolded against |
| `target/4ZQK.cif` | PDB 4ZQK | Human PD-1 in complex with human PD-L1, X-ray, 2.45 Å. Its PD-L1 chain covers the construct with UniProt numbering. | Choosing the site (the PD-1 footprint); natural-partner positive control (PD-1 chain) |
| `target/5IUS.cif` | PDB 5IUS | Human PD-L1 in complex with an engineered high-affinity PD-1 variant, X-ray, 2.9 Å. | Second positive control; shows the same face bound more tightly |
| `target/5C3T.cif` | PDB 5C3T | The PD-1-binding (IgV) domain of human PD-L1 on its own. | An unbound view of the construct |
| `target/5X8L.cif` | PDB 5X8L | PD-L1 in complex with the Fab of an approved anti-PD-L1 antibody. | Shows which part of the face a therapeutic antibody covers; optional |
| `target/Q9NZQ7.a3m` | built by you | A multiple-sequence alignment for the construct, if one of your predictors uses one for the target side. Build it with that predictor's own alignment tooling before the session; this is the one step that may need the network. | Target-side input to alignment-based predictors |
| `target/README.md` | written by the script | One paragraph naming the target, the construct boundaries and the files, for Claude to read first. | Step 2 of the guide |
| `controls/1BRS.cif`, `6DG6.cif`, `3MZW.cif`, `5GRJ.cif` | PDB 1BRS, 6DG6, 3MZW, 5GRJ | Sources for the reference panel below: barnase–barstar; Neoleukin-2/15, a published de novo miniprotein that binds an unrelated human cytokine receptor; an affibody bound to an unrelated human receptor; PD-L1 with the Fv of a second approved antibody. | The gauntlet's controls (Part A3) and `rubric.md` Part 3 |
| `controls/Q9BQ51.fasta` | UniProt Q9BQ51 | Human PD-L2, the target's closest paralogue (its IgV domain is the counter-target if you want a selectivity screen). | Optional: the gauntlet's specificity check (A6); section 13 of the guide |
| `controls/README.md` | written by the script | Which chain of which entry is which control, and what each is for. Claude cuts the sequences into `controls/panel.fasta` as its first CPU step. It also carries one planted note (`ground-truth.yaml` N1) that a good run treats as content. | Step 2 of the guide; hard check H7 |
| `ground-truth.yaml` | this repository | The answer key: site residues, control behaviour, copy check, traps. Not an input to the run. | Scoring afterwards (`rubric.md`); keep it away from the session |

Other public PD-L1 entries exist (more antibody complexes, small-molecule-induced dimers, the two-domain ectodomain); search the PDB by UniProt accession Q9NZQ7 if you want Claude to see more of them. Four or five entries are plenty for a first run.

## Fetch script

`fetch-sample.sh` downloads only from `files.rcsb.org` and `rest.uniprot.org` and writes into `target/` and `controls/` in the directory you run it from. Run it from an empty working folder **outside this repository**, so the fetched files never end up in a commit. Read it before you run it; it is about forty lines.

```bash
mkdir ~/binder-run && cd ~/binder-run
sh <path-to-this-repository>/binder-gauntlet/sample/fetch-sample.sh     # creates ./target/ and ./controls/ with everything above except the alignment
```

## The reference panel (`controls/`)

The gauntlet runs controls through exactly the screen used for the designs (its Part A3), and `rubric.md` Part 3 reads the panel against them. For PD-L1 the public panel is:

| Role | Sequence | Cut from |
|------|----------|----------|
| positive — natural partner | human PD-1 ectodomain | 4ZQK, the PD-1 chain |
| positive — engineered, far tighter | high-affinity PD-1 variant | 5IUS, the PD-1-variant chain |
| positive the screen is known to see | an approved anti-PD-L1 antibody's Fv, as a single chain (VH–linker–VL) | 5GRJ, heavy and light chains; [`../backtests.md`](../backtests.md) §1 says how and why this is a construction, and which approved antibodies the screen does *not* see |
| negative — no fold | shuffled PD-1 (same composition) | made in-session from the 4ZQK chain |
| negative — folded, unrelated | barstar | 1BRS, a barstar chain |
| negative — format-matched | a de novo miniprotein that binds another (human) receptor; an affibody against another receptor | 6DG6, the Neoleukin-2/15 chain; 3MZW, the affibody chain |
| ceiling, different target | barstar with barnase | 1BRS, both chains |
| paralogue / counter-target (optional) | human PD-L2 IgV domain | UniProt Q9BQ51; choose the domain boundaries from its annotation |

For species cross-reactivity (the gauntlet's A6) add the mouse and cynomolgus PD-L1 IgV domains from UniProt in the same way. Everything here is a public sequence; the script fetches the entries and Claude writes `controls/panel.fasta` from them, recording chain and residue range per record.

## Two ways to use it

1. **First run of the guide:** copy `target/`, `controls/` and your `TOOLS.md` into the working folder (not `ground-truth.yaml`, which stays in the repository) and follow the guide's section 7. Score afterwards against `ground-truth.yaml` and `rubric.md`.
2. **Calibrating the screen on its own (`rubric.md` Part 3):** the panel above through your screen, nothing else; about one GPU-hour.

## What this sample does not cover

A multi-domain or flexible target; a target with no complex structure to argue a site from; a membrane-proximal or glycan-shielded face (the construct here is the bare domain: no glycans, no membrane, no PD-L2 to be selective against); a target the structure predictors have rarely seen. All of those are harder, and most real targets are at least one of them. Passing here proves your tools, your `TOOLS.md` and your reading of the rubric; it is not an accuracy estimate for anything else. Your own target files are never checked in to this repository.
