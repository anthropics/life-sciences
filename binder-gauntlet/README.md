# Binder Gauntlet: design and check de novo miniprotein binders

**From a public target structure to a scored, checked design panel**

## 0. Quickstart

Run your first design panel on a public sample target. You need Claude Code on a machine with a GPU, and the design tools from section 4 installed.

**1. Install the skill.** In Claude Code, run:

```
/plugin marketplace add anthropics/life-sciences
/plugin install binder-gauntlet@life-sciences
```

**2. Make a working folder and fetch the sample target.** In a terminal, run:

```bash
git clone https://github.com/anthropics/life-sciences ~/life-sciences
mkdir ~/binder-run && cd ~/binder-run
sh ~/life-sciences/binder-gauntlet/sample/fetch-sample.sh
cp ~/life-sciences/binder-gauntlet/TOOLS.md.example TOOLS.md
```

**3. Fill in `TOOLS.md`.** Tell Claude how to call each tool on your machine. This is the only file you write.

**4. Check your setup.** Start Claude Code in `~/binder-run` and paste the setup check from section 7, step 2. Fix anything it reports.

**5. Run the design.** Paste the design brief from section 7, step 3. A first run takes about an hour on one GPU.

**6. Review the results.** Open `report.md`, `designs.csv` and `dossier.md`. Section 9 explains what good looks like and what to check by hand.

## 1. Aim

This guide takes you from a target structure to a shortlist of binder designs that is ready for design review.

You start with one target domain that has a public structure. You finish with two things:

**A panel of about 30 de novo miniprotein designs of 50 to 120 residues.** Each design comes with:

- its sequence
- a model of the designed complex
- the site it targets and the reason for choosing that site
- interface-confidence scores from at least two structure-prediction programs
- a check that the programs agree on the pose
- a structural-cluster label

**A validation dossier.** The dossier contains:

- the same screen, run on known binders and non-binders of the target
- a developability profile of every design, shown against reference proteins
- the candidates you set aside, and why
- a list of what you did not check

You can trace every number to the program that produced it. You can then take the shortlist to the people who will express and assay it.

## 2. Who it is for

This guide is for protein engineers, computational biologists and structural biologists in discovery groups. These readers already know what a binder campaign involves. They want to run the full design, screen and check loop, with the reasoning written down.

You should be able to:

- read a protein structure
- interpret ipTM, pLDDT and interface predicted-aligned-error scores
- judge the result, because a confident prediction does not prove that a design binds

This guide does not teach protein design.

You also need:

- **GPU compute that you control.** One GPU machine is enough for the first run. To scale out, use a Kubernetes namespace or a Slurm partition with GPUs.
- **Someone who can install the tools.** This person installs open-source structure-prediction and design tools on the machine or in a container image (section 4).

This guide does not suit:

- anyone who wants a sequence to order without reviewing the models
- targets with no experimental structure of the domain of interest
- antibody or peptide design, which use different tools
- any target or purpose outside your declared use case (section 11)

## 3. Surface

**Start with Claude Code and the `binder-gauntlet` skill on a GPU machine you control (Route A in section 7).** This can be a workstation, a cloud instance in your own tenant, or an interactive session on a cluster node.

Claude Code manages the loop. It chooses the site, runs the tools, reads the results, applies the checks and writes the report. The design and prediction tools run on the same machine, so you have nothing to set up between the session and the GPU.

The `binder-gauntlet` skill adds validation and developability checks that steer selection during the run.


You need only Claude Code and the open-source tools. Your compute, files and designs stay inside your own infrastructure. We developed and tested the method on this setup, so use it for your first run.

**Scale out by submitting each GPU step as a batch job to your Kubernetes or Slurm cluster (Route B).** Claude Code runs wherever `kubectl` or `sbatch` runs. It writes and submits one job per tool per batch, reads the logs and result tables, and decides what to run next. The brief and the checks stay the same.

[`cluster/README.md`](cluster/README.md) describes what the cluster needs, gives an example Job manifest and lists the rules the session follows. We have documented this route, but it has not yet been run on a customer's own cluster. Move to it after Route A has proven your `TOOLS.md`.

## 4. Prerequisites

**Access**

- Your organisation is approved under the program.
- You belong to the group that holds the grant.
- The program's model appears in the model selector on your surface.
- Your declared use case covers this work (section 5).

In Claude Code, run `/status` to confirm the organisation. Then select the program's model with `/model <id>`.

**Compute**

For the first run, you need a machine that can run Claude Code and has:

- one modern data-centre GPU with at least 40 GB of memory (smaller cards work with smaller batches, but take longer)
- about 8 to 24 CPU cores
- 100 GB of disk for model weights and intermediate files

To scale out later (Route B), you need one of the following:

- a Kubernetes namespace with a GPU node pool, where you can create Jobs, Pods and a PersistentVolumeClaim
- a Slurm partition with GPUs

You do not need cluster-admin rights. [`cluster/README.md`](cluster/README.md) gives the exact requirements.

**Software**

Install the tools before the session starts, either on the machine (Route A) or in the container image your jobs use (Route B). Install each tool from its own repository, and check its licence against your organisation's policy.

| Group | Tools |
|-------|-------|
| Backbone generators | At least two of: RFdiffusion or its successors, BindCraft, or another open de novo binder generator your group trusts |
| Sequence design | ProteinMPNN or SolubleMPNN (the LigandMPNN code base) |
| Structure prediction for screening | At least two of: AlphaFold2-multimer (ColabFold), Boltz-2, Chai-1, ESMFold2, Protenix. Claude uses two during design. Ideally, hold back one or two more for the final check (section 9). |
| Analysis | TM-align or US-align, DockQ, DSSP, freesasa or another solvent-accessibility library, Biopython, and an ipSAE script (the interface predicted-aligned-error score of Dunbrack, 2025) |

Anthropic's [September 2026 research post](https://www.anthropic.com/research/claude-uplifts-biomolecular-modeling) links to a public release of accelerated builds at [`github.com/anthropics/uplifting-biomolecular-modeling`](https://github.com/anthropics/uplifting-biomolecular-modeling). The release has one kit per tool. It covers RFdiffusion, ProteinMPNN, ColabFold, Boltz-2, Chai-1, Protenix and ESMFold2, including ESMFold2's fast, alignment-conditioned and no-alignment variants.

The kit code uses the Apache-2.0 licence. Each upstream tool keeps its own licence, which you will find in the kit's folder. The release is published as-is and is not maintained. Check each upstream licence as you would for the tool itself, and record the builds you use in your run notes.

Write a one-page `TOOLS.md` from [`TOOLS.md.example`](TOOLS.md.example). It should say where the GPU work runs, how to call each tool, where the weights are and how much GPU memory each tool needs. Claude reads this file first and writes every tool call from it.

**The plugin and the skill**

1. In Claude Code, run `/plugin marketplace add anthropics/life-sciences`.
2. Run `/plugin install binder-gauntlet@life-sciences`.
3. Ask "which skills do you have?" and confirm that the list includes `binder-gauntlet`.
4. Run the skill's self-test once: `python3 ~/life-sciences/binder-gauntlet/scripts/sequence_profile.py --self-test`.

**Network**

The run needs no network access from the tools. Fetch the target files (section 6) and model weights in advance. On a cluster, bake the weights into the image. The session itself still sends what it reads to the model service (section 10).

**Time**

| Task | Time |
|------|------|
| Tool installation (once, by your compute administrator) | Half a day to a day |
| First run on the public sample target, on one local GPU | About an hour of session time, or nearer two on a cluster, where each step waits to schedule |
| Controls and held-out check (within the session time above) | About one GPU-hour in total |
| Your own reading, scoring and hand checks | Two to three hours |
| A serious run on your own target | A working day of session time and a day of review |

## 5. Declared use case

Paste and adapt the statement in [`DECLARED-USE-CASE.md`](DECLARED-USE-CASE.md). In short:

> Computational design and in-silico screening of de novo miniprotein binders against named protein targets in our discovery portfolio, by our protein-engineering and computational-biology staff, using open-source design and structure-prediction tools driven by the program's model on compute we control; public or our own target structures in, a scored and clustered design panel with a validation and developability dossier and a written rationale out; every panel reviewed by a protein scientist before anything is ordered or expressed; no toxins, no pathogen proteins, no patient data.

Make the statement specific to your actual targets and team. A specific statement makes monitoring simpler.

## 6. What you bring

| Input | Format | Notes | Required |
|-------|--------|-------|----------|
| Target domain sequence and boundaries | UniProt accession and residue range; FASTA | Use the exact construct you would express for an assay. Every score uses this chain, so decide it first. | Yes |
| Experimental structures containing the domain | mmCIF or PDB files from the Protein Data Bank | Include complexes with natural partners or existing binders if they exist. They show which surfaces are functional. They also supply the **positive controls** that the gauntlet runs through the same screen as your designs. | Yes |
| Reference panel for the checks | FASTA with a few sequences, each with its public source | Positives: the natural partner or a known binder, cut from a deposited complex. Negatives: a shuffled positive, one or two unrelated well-folded proteins of similar size, and one or two format-matched binders of *other* targets. Optional: the target's closest human paralogue and the mouse or cynomolgus orthologue, for specificity and cross-reactivity checks. [`sample/README.md`](sample/README.md) lists the public PD-L1 sequences. [`backtests.md`](backtests.md) shows how they perform in the screen. | Recommended |
| Counter-target | Same as the target | Provide a close homologue that the binder should *not* bind, if selectivity matters. Claude refolds candidates against both (section 13). | Optional |
| Multiple-sequence alignment for the target | a3m | Some predictors score the target side better with an alignment. Build it in advance with the predictor's own tools. The binder side is always a designed sequence with no alignment. | Recommended |
| `TOOLS.md` | Markdown, one page | States where the GPU work runs and how to call each tool (section 4). | Yes |
| What you already know | A few sentences in the prompt | Tell Claude which surface matters biologically, which epitopes your existing molecules cover, any length or composition limits from your expression system, and the intended format (bare binder, Fc fusion or tracer). | Optional |

**On your first run, bring nothing of your own.** Use the public sample target in [`sample/README.md`](sample/README.md): the IgV domain of human PD-L1 (UniProt Q9NZQ7, residues 18–132). The sample includes a small public reference panel and four public Protein Data Bank entries:

- the domain alone
- its complex with PD-1
- its complex with an engineered high-affinity PD-1
- its complex with a therapeutic antibody fragment

PD-L1 is a well-studied, rigid, single-domain target with a known functional face. This makes it a fair test of your setup and an easy case for the method. Your own targets will probably be harder.

## 7. Step by step

The prompts below are plain text. Use them unchanged on your first run. On later runs, change the target lines, the constraints and the budget, but keep the structure. The design brief is the same on both doors. Only the setup differs.

### Route A: Claude Code on a GPU machine you control

**Step 1. Install and select (5 minutes).** In Claude Code, on the GPU machine, run:

- `/plugin marketplace add anthropics/life-sciences`
- `/plugin install binder-gauntlet@life-sciences`
- `/status`, and confirm that the organisation is correct
- `/model <program model id>`

**Step 2. Prepare the working folder and check the tools (10 to 20 minutes the first time).** Most of this time goes to correcting `TOOLS.md`.

1. Create a working folder **outside this repository** and start Claude Code in it.
2. Add your `TOOLS.md`, plus the `target/` and `controls/` folders that the fetch script in [`sample/README.md`](sample/README.md) builds.
3. Paste this setup check:

   ```
   Before we start, check my setup. Confirm this folder holds target/ (with README.md, construct.fasta
   and at least one structure file), controls/ and a filled-in TOOLS.md, and list anything missing.
   Then run nvidia-smi and the smoke test TOOLS.md gives for each tool, one at a time. Show me a short
   table of the GPU and each tool (OK, or failed with the first error line), then stop.
   ```

4. Wait until every tool in `TOOLS.md` has responded on this machine before you continue.

**Step 3. Give the design brief (about an hour of session time on one local GPU).** Paste this brief:

```
I would like a panel of de novo miniprotein binder candidates against the target described in
target/README.md (human PD-L1, IgV domain, UniProt Q9NZQ7 residues 18-132; the structures in
target/ are public Protein Data Bank entries; controls/ holds a small public reference panel with
a README saying what each sequence is). TOOLS.md says where the GPU work runs and how each
installed tool is called; batch candidates so each tool run loads a model once, and keep every
script you write under jobs/ with its log. Plan for roughly two hours of my wall-clock time on
one GPU at a time (tell me if I may use more); keep your own eye on the time, and leave the
last fifteen minutes for tidying up and final checks.

Please work like this, and tell me briefly what you are doing as you go:
1. Look at the structures and choose a binding site. Say which surface you chose, which
   residues, and why a binder there would do something useful.
2. Run the reference panel in controls/ through exactly the screen you will use for candidates,
   first, and write the numbers to controls/reference-points.md, so every later number has
   something to be read against.
3. Generate candidate backbones against the site with at least two different generators, and
   design sequences onto them. Aim for a few hundred candidates of 50-120 residues.
4. Refold every candidate in complex with the target using two structure predictors. Keep a
   candidate only if both are confident about the interface AND their refolded pose agrees with
   your design model. Record the numbers per candidate.
5. Use the binder-gauntlet skill at every selection round, not just at the end: run its cheap
   checks on the current shortlist and let them change what you keep, drop or redesign; keep a
   running list of what you set aside and why.
6. Cluster the survivors by backbone shape and pick about thirty that cover as many clusters as
   the scores allow. Check that no sequence is a close copy of the target itself or of anything
   in the structures or controls I gave you.
7. Write designs.csv (one row per design: id, sequence, length, generator, site, each
   predictor's interface confidence, pose agreement, cluster, one line of notes), one model file
   per design under models/, report.md (the site and why, what you ran, how you chose and
   steered the thirty, what the numbers do and do not mean, what a protein scientist should
   verify before ordering anything), the gauntlet's dossier.md with its per-design CSV, and
   the scripts you wrote under jobs/.
Save a first complete version of designs.csv, models/ and report.md early and keep them
current, so that if anything stops the session I still have a usable panel.
```

Claude reads `TOOLS.md` and the structures, and reports the chosen site within a few minutes. It then runs the controls through the screen and starts the first generator batches. You do not need to watch.

**Step 4. Read along, or come back later.** If you read the transcript, look for the following:

- Did the controls go through the screen before the candidates?
- Did Claude notice disagreement between predictors, and what did it do about it?
- Did the gauntlet's checks change a selection? For example, Claude might drop a candidate for a deamidation site in the interface, or redesign one with a hydrophobic free face. Or did Claude only write up the checks at the end?
- Did Claude write the files early?

If Claude asks a question, answer it. Otherwise, let it work.

**Step 5. Ask for the self-check when Claude says it has finished (5 minutes).** Paste:

```
Before I read the panel: re-open designs.csv, report.md and dossier.md and check them against
what actually ran. For each design confirm the model file exists and the sequence in it matches
the row. List any number in report.md or dossier.md that you did not compute in this session,
and anything the dossier's "not checked" section is missing. Then tell me the three designs you
are least sure about and why.
```

**Step 6. Run the held-out check (about half a GPU-hour).** Ask Claude to refold the final panel and the controls with the predictor you kept out of the session. Ask it to add the result to section 1 of the dossier as the independent check. Then drop any design that this check disputes.

**Step 7. Score the run (two to three hours of your time).**

1. Score the run with [`rubric.md`](rubric.md) Parts 1 and 2.
2. On the sample target, also score it against [`sample/ground-truth.yaml`](sample/ground-truth.yaml), which covers site residues, control behaviour, the copy check and the planted note.
3. Read the dossier against Part 3 of the rubric and [`backtests.md`](backtests.md) §2.
4. Complete the hand checks in section 9 before the panel goes anywhere.

### Route B: Claude Code with GPU steps as batch jobs on your cluster

Complete Route A once first. Then run Claude Code where `kubectl` (with a context for your cluster) or `sbatch` (on a login node) runs. Make sure the image or environment described in [`cluster/README.md`](cluster/README.md) is in place.

1. **Install and select**, as in Route A step 1.
2. **Prepare the working folder (15 to 30 minutes the first time).** Follow Route A step 2, and fill in the Kubernetes or Slurm block of `TOOLS.md`. On Slurm, use the shared filesystem. Ask Claude to read that block, set up the work volume (on Kubernetes), and submit one short test job that prints `nvidia-smi` and each tool's version, then show you the result and stop. Wait until the job has run on a GPU node and every tool has responded. Most of this step goes to finding out what your cluster calls its GPUs.
3. **Send the Route A brief unchanged.** `TOOLS.md` already says where the GPU work goes. Add "follow the rules in binder-gauntlet/cluster/README.md". Continue from Route A step 4, and expect the run to take nearer two hours, because each step waits to schedule.
4. **Clean up.** Ask Claude to remove what it created for this run. It removes only objects that carry this run's labels, by name. It keeps the volume until you tell it otherwise.

## 8. What you get

You get the following files in the working folder. On Route B, Claude copies the tables and selected models back from the work volume as the run progresses.

| File | Contents |
|------|----------|
| `designs.csv` | One row per design: `id, sequence, length, generator, site`, a `<predictor>_interface_confidence` column for each predictor, `pose_agreement` (DockQ between the design model and each predictor's refold, or their mean), `cluster` and `notes` |
| `models/` | One PDB or mmCIF file per design, showing the designed complex as it came out of the pipeline, with binder and target as separate chains. These files are design models. They are not later predictions or experimental structures. |
| `report.md` | The chosen site and the reasoning; the tools and settings used; the number of candidates generated, screened and kept at each stage; where the gauntlet's checks changed a decision; how the final panel was selected; the numbers and what they measure; and a list of items to verify by hand |
| `dossier.md` and `dossier.csv` | The `binder-gauntlet` skill's validation dossier, plus one row per design and control, keyed by sequence so it joins to `designs.csv`. The format is in the skill's `references/dossier-template.md`. |
| `controls/reference-points.md` | The control sequences' scores through your screen. Part 3 of [`rubric.md`](rubric.md) asks for this table. |
| `jobs/` | Every script Claude wrote for a GPU step, with its log (and on Kubernetes, the manifest it applied). Together with the session transcript, these files record what ran. Keep them with the outputs. |

The dossier has nine sections:

1. Binding evidence beyond the objective, with the controls' numbers beside the designs'
2. Function or blocker activity
3. Specificity and species cross-reactivity
4. Robustness
5. Native context
6. Developability per design
7. Interface metrics
8. Candidates set aside, and why
9. What was not checked, and why

## 9. What good looks like

This section summarises [`rubric.md`](rubric.md).

**The structures support the site choice.** The report names residues that exist in the deposited entry and says what binding there would do. "A hydrophobic patch" is not a reason. "The face PD-1 engages in entry 4ZQK, so a binder there competes with PD-1" is a reason.

**Every kept design passed at least two predictors and a pose-agreement check.** The report shows the numbers per design and per predictor. A panel screened on one predictor inherits that predictor's blind spots.

**The controls went through the identical screen first.** Their numbers sit beside the panel's numbers in the dossier. If the screen cannot separate the target's known partner from a shuffled sequence, the panel's numbers tell you little.

**The checks steered the run.** The dossier's "set aside" section has entries. Claude dropped or redesigned some candidates mid-run because of an interface liability, a hydrophobic free face, a self-association hit or a failed monomer refold, and the report says so. If Claude wrote the dossier only after the final pick, the dossier says that too.

**The panel is structurally varied.** The report says how Claude measured variety. Thirty sequences on one backbone are one bet, not thirty.

**No design is a copy.** No design shares long stretches of sequence with the target, its natural partners or any binder in the structures or controls you supplied. The report says how Claude checked this. A wider novelty search is a hand check (see below).

**The report describes each number accurately.** Interface scores measure predictor confidence in a modelled interface. Developability numbers are descriptors shown against reference proteins, with no pass or fail where the literature gives no threshold. The report does not call any design a binder or rank designs by "affinity". The last section of the dossier lists what was not assessed.

**The back-tests set your expectations.** See [`backtests.md`](backtests.md).

On PD-L1, a three-program screen gave these results:

- It kept all thirteen public negatives below 0.25.
- It found natural and engineered PD-1.
- It **missed the binding domains of three of the four marketed anti-PD-L1 antibodies tested**.
- It could not tell a micromolar partner from a picomolar one.

Expect designed helical miniproteins to score above real binders on this screen, because of how they are built.

For developability, §2 gives a public reference panel of small-protein drugs with computed descriptor ranges. It also names three things to catch in round one for helical de novo designs:

- no tryptophan or tyrosine, which rules out an A280 measurement
- deamidation-prone Asn and Asp motifs in the interface
- alanine-rich stretches that read as hydrophobic on sequence scales

Read your own dossier against these results.

**Check these by hand before anyone orders or expresses a design, or shows it outside the team:**

1. Open several design models next to the target structure. Confirm that the binder sits on the claimed site with a plausible interface.
2. Search every sequence against public databases and your own collections for close matches to known proteins and existing binders of this target.
3. Read the dossier's developability table against your expression system's rules.
4. Read the "not checked" section against your own diligence list.
5. Confirm that the team applied the held-out predictor's verdict.

The designs that survive are candidates for expression and a binding assay. A qualified protein scientist decides what goes to the bench. The assay is the first point at which anyone can say that a design binds.

## 10. Data handling

**What we send.** The session sends the following to the model on your chosen surface, under your organisation's program terms:

- your prompts
- the files Claude opens or writes, including target files, `TOOLS.md`, job scripts, result tables, the report and the dossier
- any job output that Claude reads

The design and prediction tools run on your own machine or cluster. They make no network calls once their weights are in place. On a cluster, the jobs need no egress. The work volume leaves the cluster only as the text and tables that the session reads back.

**Where it is processed and how long it is kept.** Processing takes place on Anthropic's first-party service, under the program's published terms.

- Anthropic retains program usage for 30 days on every program surface.
- Zero data retention is not available. Existing zero-retention or custom-retention agreements do not apply to usage under a grant.
- Automated systems screen usage. People at Anthropic review flagged activity only, in rare circumstances.
- Anthropic never uses program usage for model training or product improvement.
- Anthropic monitors usage against your declared scope. If something looks out of scope, Anthropic notifies your admin and gives you one week to explain, narrow the scope or dispute it.

**What stays with you.** Claude writes every output file to your working folder or work volume, and writes nothing back to this repository. The machine or image, the weights, the volume and the jobs all remain yours. On a cluster, the session creates objects only in the namespace you name. It labels everything it creates, so that it can remove only those objects ([`cluster/README.md`](cluster/README.md)).

**What is sensitive.** On the first run, nothing is sensitive: the target is public and the designs are new sequences that belong to you. On your own targets, the target identity, the site and the designed sequences may all be confidential to your programme. Keep the work volume, the image registry and the outputs where such material is allowed. Remember that the transcript and the job logs contain this material too.

**What not to send.** Do not send:

- patient-derived sequences or data
- protected health information (PHI)
- sequences or structures you do not have the right to process
- anything outside your declared use case

During the beta, no Business Associate Agreement (BAA) covers the program. PHI workloads are therefore out of scope, and an organisation set up under a BAA cannot receive a grant.

## 11. Access level

De novo binder design against a human therapeutic target with open-source tools is mainstream discovery computational biology. It runs on the program's **Standard Use** grant. The closest published categories are:

- "Drug and countermeasure discovery and preclinical support"
- "Bioinformatics, computational biology, and research data infrastructure"

The program's launch post names basic science and R&D among the areas Standard Use serves. The use case your organisation declared and had approved decides whether your own work is in scope. Ask your program contact if you are unsure which category to use.

Run this work inside your program organisation with the program's model selected. Do not use a personal account or a non-program model. Outside a grant, the standard safeguards on biological research may decline protein-design sessions like this one.

This guide does not need the High-risk Use add-on. Toxins, pathogen proteins and anything your institution handles under enhanced containment are a different project. Ask your program contact before you declare such work, and do not run it under this guide's statement.

If the model declines a step that you believe your declared use case covers, do not rephrase the request to get around it. Talk to your program contact instead.

## 12. Troubleshooting

| Symptom | Likely cause | What to do |
|---------|-------------|-----------|
| A generator or predictor crashes on start, or Claude cannot find it | The tool is not on the path, the weights are missing, or an expected scratch-directory or environment variable is unset | Update `TOOLS.md` with the exact working command and any environment variables. Make the tool check in step 2 run each tool's smallest example before the real session. Claude can usually work around failures that reach a log. It cannot work around silent ones. |
| The GPU runs out of memory, or everything slows down | Several tools are loaded at once, or batch sizes suit a larger card | Add the GPU's memory to `TOOLS.md`. Run one tool per run and batch inside the run. On small cards, drop the largest predictor from the in-session screen and use it only for the held-out check. |
| On a cluster, jobs stay **Pending** | No GPU node is free, your GPU pool needs a node selector or toleration, or your cluster names the GPU resource differently | Read the pending pod's Events lines, or ask Claude to describe the pod and read them to you. Fix the Kubernetes block of `TOOLS.md`. Ask your platform team what a GPU job in this namespace must carry. [`cluster/README.md`](cluster/README.md) covers more cases, including image pull errors, volume access modes and OOMKilled. |
| The skill never ran, or the dossier appeared only at the very end | The request did not name the skill, or the session treated it as a reporting template | Name the skill in the first message ("use the binder-gauntlet skill") and keep point 5 of the brief. A dossier written after the run still has value, and the skill says to label it that way. |
| The two predictors disagree on most candidates, or the panel is 30 near-identical bundles | Single-sequence and alignment-based predictors often disagree, or selection ranked on score alone | Keep only candidates that both predictors accept in the design pose, and report the disagreement rate. Ask for clustering by backbone shape with a cap per cluster, and ask Claude to state the number of distinct clusters. The held-out check in step 6 decides. |
| The model declines the design task or part of it | You are not on the program's model in the program organisation, the grant has not reached your group, or the work is outside your declaration | Do not rephrase to get around it. Check the surface, organisation, model and grant (section 11). Check that your declared use case names this kind of work and these targets. |

## 13. What's next

**Your own target.** Use the same brief with your construct, structures, reference panel and constraints. Run the controls first. If the screen cannot separate your target's known partner from a shuffled sequence, it will not tell you much about designs. If your known binders are antibodies, expect the screen to miss some of them ([`backtests.md`](backtests.md) §1).

**Selectivity and species cross-reactivity.** Add a counter-target (a close homologue the binder must avoid) and the orthologues your animal models need. Ask Claude to refold every shortlisted candidate against each one. The results go in section 3 of the dossier. Each extra chain adds roughly one screen's worth of Jobs.

**More GPUs and longer sessions.** The steps are independent, so several GPUs can run generators and predictors side by side, on one machine or as cluster jobs (Route B). You only need to change the budget sentence in the brief. A day of session time with interface redesign evens out the scores. Only experiments can show whether this produces better binders, so the checks matter even more as the search gets stronger.

**Physics-based and licence-gated checks.** Rosetta or FoldX energies and dedicated aggregation, solubility and T-cell-epitope predictors stay under "not checked" in the dossier until you add them to your image and `TOOLS.md`. The skill then uses them.

**To the bench.** Give the surviving panel, the report, the dossier, the models and the controls to the team that runs your expression and binding assays. Include the binding-dead mutant that the dossier names for each lead, as the assay's specificity control. Use the natural partner as a positive control. A yes or no binding readout is the first result that means anything.

## License

Apache 2.0. See [LICENSE.txt](LICENSE.txt).
