# Step catalogue — what runs, in what shape

The design-and-screen loop in the [guide](README.md) breaks into the steps below. On a GPU machine (Route A) each GPU step is one tool run in the working folder; on a cluster (Route B) it is one batch job as [`cluster/README.md`](cluster/README.md) describes; the CPU steps run wherever Claude Code runs. Commands are written generically on purpose: the exact invocation of each tool on *your* machine or in *your* image is what your `TOOLS.md` records (start from [`TOOLS.md.example`](TOOLS.md.example)); nothing here is a real flag of any tool. Paths are relative to the working folder (`/work` inside a cluster job). Keep every script you write for a step under `jobs/`, with its log beside it.

Times are orders of magnitude for a first run — about thirty designs chosen from a few hundred candidates of 50–120 residues against a single domain, on one data-centre GPU with 40 GB of memory or more — and depend on the tool, the build, the GPU and the lengths. Time your own first small batch and plan from that, not from this page.

**Batching rule.** One tool run (or job) per tool per batch, never one per candidate: loading a model takes longer than scoring one more complex, and on a cluster every job also pays scheduling and image start-up. Batches of 50–200 candidates are a good unit. Two light steps that always follow each other (backbone generation then sequence design) can share one run as `step-one && step-two`.

## 0. Target preparation — CPU

Inputs: the target FASTA, the public structures, the alignment if you built one (guide section 6). Claude Code trims the construct, picks the site and writes `target/README.md` with the residues chosen and why. On a cluster, `target/` and `TOOLS.md` are then copied onto the work volume. Output: `target/`.

## 1. Backbone generation — GPU

Inputs: `target/`, the site residues, a length range. Command shape: the generator's own command line as `TOOLS.md` gives it, writing into `backbones/<batch>/`. Use at least two different generators across the campaign (for example RFdiffusion or a successor, BindCraft, PXDesign, BoltzGen, Proteina-Complexa, the ESM cookbook's gradient-guided design through ESMFold2 or another open generator you have installed; some of these design the sequence as well, in which case step 2 is a redesign or is skipped for that batch). A generator that optimises each design against a structure predictor makes that predictor part of your objective: for example, BindCraft does this with AlphaFold2-multimer, the ESM cookbook recipe with ESMFold2 and mosaic with Boltz-2. Do not count that predictor as an independent check on that generator's designs, and hold back a different one for step 6. Output: one model per candidate backbone in complex with the target. Typical: minutes to tens of minutes per batch of 50–200; memory well inside a 40 GB card at these sizes.

## 2. Sequence design — GPU, light

Inputs: `backbones/<batch>/`. Command shape: ProteinMPNN or SolubleMPNN (both from the LigandMPNN code base) over the whole batch folder, target chain fixed, a few sequences per backbone, writing `sequences/<batch>.fasta` plus a table mapping sequence id to backbone. Typical: seconds to a few minutes per batch. Cheap enough to share a run with step 1.

## 3. Complex re-prediction — GPU, the bulk of the time

Inputs: `sequences/<batch>.fasta`, the target sequence (and its alignment for predictors that use one). Command shape: one predictor over one batch, each candidate as a two-chain complex with the target, *N* seeds, writing models and confidence files under `predictions/<predictor>/<batch>/`. Run two or three predictors that see the problem differently — whichever of ESMFold2, Protenix, AlphaFold2-multimer (through ColabFold), Boltz-2 and Chai-1 you have installed — and keep one more back for step 6. **One run per predictor per batch.** Typical: of the order of 10–60 seconds per complex per seed, so this step dominates; screen cheaply first (one fast predictor, one seed, the whole batch), then spend seeds and the slower predictors only on what survives. Memory: the alignment-based and diffusion-based predictors are the ones that need the larger cards and more host memory.

## 4. Scoring, pose agreement and clustering — CPU

Inputs: `predictions/`, `backbones/`. Command shape: an ipSAE script over each prediction's aligned-error output (the interface score of Dunbrack, 2025), ipTM and interface pLDDT from the predictors' own confidence files, DockQ or an interface RMSD between the design model and each refold (pose agreement), and TM-align or US-align all-against-all on the binders for structural clusters. On a cluster this runs in the utility pod or as a CPU-only job, or locally after copying the confidence files back. Output: `scores/<batch>.csv`, one row per candidate per predictor, merged into `designs.csv`. Typical: minutes.

## 5. In-loop checks from the binder-gauntlet skill — CPU, plus more of step 3

The [`binder-gauntlet`](SKILL.md) skill's cheap checks (sequence composition, similarity to the target and its known partners, liabilities, and the like) run where Claude Code runs, on `designs.csv`, every selection round, and cost no GPU time. The structural controls it asks for — the binder folded alone, a shuffled-sequence control, the natural partner as a positive reference — are further step-3 runs over small FASTA files; batch them together rather than one run each.

## 6. Final pass with the held-out predictor — GPU

Inputs: the thirty or so selected designs as one FASTA. Command shape: the predictor you did *not* use in steps 3–5, more seeds than the screen used, writing `predictions/final/`; then step 4 once more over that folder. One run. Typical: well under an hour at this size. Its numbers go into the report next to the screening numbers, with the disagreement stated, as the guide's section 9 asks.

## Parallelism

On one GPU, run one GPU step at a time and do CPU work (reading earlier results, the skill's checks, the report) while it runs. On a machine with several GPUs, or on a cluster, steps are independent: start several in the background — different batches, or different predictors over the same batch — up to the number of GPUs you have or your quota allows. On a cluster, more than that only queues jobs as `Pending` while their deadlines tick.
