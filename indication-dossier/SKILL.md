---
name: indication-dossier
description: >
  Write a cited dossier on one therapeutic indication: its patients,
  epidemiology, biology, standard of care, regulatory path and key trials. Use
  when assessing a disease area for drug development or trial design, or when
  asked to "write an indication dossier", "get me up to speed on [indication]"
  or "characterize this patient population". Input is an indication name.
---

# Indication Dossier

Produces a structured research dossier on a single indication, framed as a
patient population: who they are, what drives their disease, how they're
treated today, and how clinical trials can be designed to help them. Runs as
five phases that write resumable waypoint files; after a brief identity check
at the end of Phase 1, the remaining phases run straight through.

## Framing

**Think of an indication as a patient population.** Frame everything
from the patient perspective: "Who are these patients?"; "How are these
patients identified and managed?"; population nesting: "all patients in
{child} are patients in {parent}".

Some indications don't map to ICD codes or standard disease definitions:
"immunosenescence" is a biological state, not a billable diagnosis; "ageing"
is not an FDA-accepted indication; "GLP-1 induced sarcopenia" is an
iatrogenic population. Note these distinctions explicitly. They matter for
regulatory path and trial design.

## Inputs

- **`indication`** (required) — indication name (e.g., "sarcopenia",
  "idiopathic pulmonary fibrosis").
- **`additional_context`** (optional) — areas to focus on, parent
  indication, or other framing.
- **`workdir`** (optional) — where the waypoints and the final report go.
  Defaults to `indication-dossier-<slug>/`.

## Tools this skill expects

| Purpose | Tool |
|---|---|
| ClinicalTrials.gov | A clinical trials MCP server or connector, if available |
| Literature | A PubMed MCP server or connector, if available |
| Web | Web search and web fetch: FDA and EMA guidance and approval pages, treatment guidelines (NCCN, AASLD, specialty societies), CDC and WHO epidemiology pages |
| Documents | Web fetch reads a remote page or PDF as text; see the research standards for PDFs whose figures carry the data |
| Parallel research | If a subagent tool is available, use it for Phases 2-4: one subagent per source family in a phase, each given its part of the phase's reference file, the research standards and `meta.json`. Tell each one to return findings rather than write waypoint files, and to return with them every URL it tried and how far it got (read in full, excerpt only, failed) and every caveat about a source. Running the same searches yourself, in turn, is the fallback for when no subagent tool exists; it covers less ground |

If a clinical trials or PubMed connector isn't available, say so and fall back to web search against the public site (clinicaltrials.gov, pubmed.ncbi.nlm.nih.gov), recording what the fallback loses.

## Output layout

```text
<workdir>/
└── waypoints/
    ├── progress.json                 # where the run has got to
    ├── meta.json                     # phase 1
    ├── epidemiology.json             # phase 2
    ├── biology_soc.json              # phase 3
    ├── regulatory_trials.json        # phase 4
    ├── sources_evaluated.json
    ├── research_output.json          # phase 5 — structured output
    └── indication_dossier_report.md  # phase 5 — the deliverable
```

Schemas for every waypoint file are in `references/waypoint-schemas.md`.
Waypoints are the resumable state; see "Resuming" below.

## Before starting

Read `references/00-research-standards.md`. It governs sourcing, the
anti-fabrication rules, and how to handle untrusted fetched content for every
phase. Then create `<workdir>/waypoints/` and `progress.json`, and update `progress.json` after each phase.

## Workflow

The dossier is built in five phases. Only the main session writes waypoint
files. The research for Phases 2-4 may run at the same time, but write the
waypoints in order (2, 3, 4), because biomarkers (Phase 3) inform endpoints
(Phase 4). Each time a waypoint is written: copy every URL that you or a subagent tried
in that phase into `sources_evaluated.json` with its result, failures included (the log
is the full list, not a selection); set `current_phase` in
`progress.json` to the next phase; and give a ≤200-word summary of what you
found and what's uncertain. Then go straight on.

Apart from asking before a download (see the research standards), the one
stop is after Phase 1. When `meta.json` is written, show the resolved
identity (name and aliases, codes, parent, what is in and out of scope; this
may run past 200 words) and ask the user to **proceed**, **revise the
identity** or **stop**. Ask even when the name is unambiguous: this is also
where the user narrows the population. Skip the question only if the request
said to run without stopping, or already named both the parent indication and
what is in and out of scope, or if the session cannot receive a reply (a
scheduled or headless run). In those cases state the interpretation you are
proceeding on and continue.

### Phase 1 — Meta initialization

Read `references/01-meta-initialization.md`. Resolve the indication identity:
clinical definition, ICD codes, aliases, parent indication, and whether it's
a recognized diagnostic entity. Run a quick registry landscape scan. Stand up
`waypoints/meta.json`.

### Phase 2 — Epidemiology research

Read `references/02-epidemiology-research.md`. Characterize the population:
diagnostic criteria, prevalence and incidence, demographics and risk factors,
natural history. Search the literature and the web side by side. Write `waypoints/epidemiology.json`.

### Phase 3 — Biology & standard-of-care research

Read `references/03-biology-soc-research.md`. Establish pathophysiology,
biomarkers, approved therapies, treatment guidelines, and unmet need. Three source families side by side:
literature for biology, web for guidelines, FDA pages for
approvals. Write `waypoints/biology_soc.json`.

### Phase 4 — Regulatory & trials research

Read `references/04-regulatory-trials-research.md`. Establish FDA/EMA
accepted endpoints, regulatory precedents, typical trial design parameters,
landmark trials, and notable failures. Three source families side by side: FDA pages for
guidance and approvals, the registry for trial patterns, literature for trial-history
reviews. Write `waypoints/regulatory_trials.json`.

### Phase 5 — Synthesis

Read `references/05-synthesis.md` and `references/06-writing-style.md`. Read
all four consolidated waypoint files. Write
`waypoints/indication_dossier_report.md` — narrative sections in the order
the synthesis reference specifies, with inline citations per the style guide
— and `waypoints/research_output.json`. Before writing, list the single values
still missing from waypoint fields (an approval year, an NCT ID, a primary
endpoint, a decision date) and try one targeted look-up for each. Do not open
new research threads in this phase: anything broader than a single value, name
in section 6 rather than filling it.

## Resuming

If invoked with a `workdir` that already contains waypoints: list which phases
are complete (waypoint file exists and is non-empty), show the meta summary,
and ask the user which phase to run next. If the session cannot receive a reply,
resume from the first incomplete phase and say so. Never overwrite an existing waypoint
without confirmation.
