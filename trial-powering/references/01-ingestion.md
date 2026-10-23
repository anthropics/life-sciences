# Phase 1: Ingest design parameters

Goal: produce `design_params.json` with everything `power_calc.py` needs,
each field tagged with where it came from.

## Target schema

```json
{
  "nct_id": "NCT...",                     // if known
  "asset": "...", "indication": "...", "phase": "...",
  "endpoint_type": "binary|continuous|tte|count|ordinal|other",
  "primary_endpoint": "text description",
  "archetype": {                          // see references/08-endpoint-archetypes.md
    "id": "A-J or other",
    "regulatory_convention": "one line: what this indication requires (e.g. cognitive + functional co-primary)",
    "convention_source": "URL of FDA disease guidance or recent approval precedent",
    "comparator_behavior": "placebo|vehicle|sham|active|natural-history, with one-line precedent response note",
    "structural_correction": "none|co-primary-joint-power|replication-set|paired-units|dilution"
  },
  "test": "superiority|non-inferiority|equivalence",
  "sided": 1 | 2,
  "alpha": 0.05,
  "stated_power": 0.80,
  "allocation": "1:1" | "2:1" | ...,
  "n_per_arm": [..],  "n_total": ...,
  "dropout_assumed": 0.10,                // fraction
  "interim_analyses": {"count": 1, "spending": "OBF|Pocock|none", "info_fractions": [0.5]},

  // endpoint-type specific (fill what applies):
  "binary": {"p_control": ..., "p_treat": ..., "delta_abs": ...},
  "continuous": {"mean_control": ..., "mean_treat": ..., "delta": ..., "sd": ..., "ancova_r2": null},
  "tte": {"hr": ..., "median_control": ..., "median_treat": ..., "events_required": ..., "accrual_months": ..., "min_followup_months": ...},
  "count": {"rate_control": ..., "rate_treat": ..., "rate_ratio": ..., "dispersion_k": ..., "exposure_years": ...},
  "ni_margin": null,

  "param_sources": {"<field>": {"source_type": "ctgov|deck|paper|sap|user", "locator": "slide 14 / table 2 / NCT field name", "url": "..."}}
}
```

Leave a field `null` and add a note in `param_sources` if it could not be
found. Do not guess load-bearing numbers.

## Ingestion path A: ClinicalTrials.gov

`clinical_trials.get_trial_details(nct_id=...)` gives you arms, allocation,
masking, target enrollment, primary outcome text, and sometimes a stat plan
synopsis under `Detailed Description` or as an attached protocol/SAP document
(check `documents` / `largeDocs`).

CT.gov almost always has `n_total` and the endpoint description. It rarely has
the assumed effect size, control arm assumption, SD, or alpha spend. Those
come from path B/C or get reverse-engineered in Phase 2.

## Deck sweep (runs regardless of ingestion path)

Even when CT.gov or a paper supplies the basics, the sponsor's presentations
are often the highest-yield public source for the statistical design.

Procedure:

1. Find the sponsor's "Events & Presentations" (or "Presentations") page.
2. Retrieve, at minimum: (a) the most recent corporate presentation, and
   (b) any data-event deck covering this trial (conference deck,
   R&D day, topline-results call deck). Prefer the data-event deck; it
   usually reproduces the medical presentation's statistical design slide.
3. Save retrievable PDFs to `<workdir>/decks/`, following the download rules
   in `00-research-standards.md`. If a deck can't be fetched, check EDGAR (the
   same deck or press release is often furnished as an 8-K exhibit), then
   record the status.
4. Repeat for every precedent trial used in the Phase 3 audit: the precedent
   sponsor's deck covering that trial's design or readout. Acquired or
   private sponsors: check the acquirer's archive and conference archives
   before recording `not_found`.
5. Write `<workdir>/deck_manifest.json` before `design_params.json`:

```json
{
  "decks": [
    {
      "company": "...",
      "role": "subject_sponsor" | "precedent_sponsor",
      "trial": "NCT... or name",
      "title": "Conference presentation",
      "date": "YYYY-MM",
      "url": "...",
      "status": "fetched | blocked | not_found",
      "local_path": "decks/....pdf or null",
      "extracted": "one line on what it yielded (e.g. 'stat design slide: alpha split, target HRs, planned events')"
    }
  ]
}
```

A memo that ships without this manifest is incomplete. Decks recorded as
`blocked` or `not_found` are acceptable outcomes; an empty manifest is not.

## Ingestion path B: Company deck / corporate presentation / SAP PDF

Use the `Read` tool on the PDF. Per `00-research-standards.md`, WebFetch is
text-only; for slide decks the powering parameters are usually on a single
"trial design" slide as a graphic, so visual reading via `Read` is required.

What to look for, by slide archetype:

- **Trial design schematic** (the swimlane / CONSORT-style box diagram): N per
  arm, randomization ratio, stratification factors, primary endpoint label.
- **Statistical considerations / powering footnote**: tiny text at the bottom
  of the design slide. This is where sponsors disclose "N=XXX provides 90%
  power to detect HR=0.YY assuming median control of ZZ months, two-sided
  alpha=0.05". Capture it verbatim as the `quote` for `param_sources`.
- **Prior data slide** (Phase 1b/2 waterfall, KM, forest): the effect size the
  sponsor is anchoring on. Record it; Phase 3 will compare it to precedent.
- **SAP / protocol** (if you have the actual document): section typically
  titled "Sample Size Justification" or "Determination of Sample Size".

If there are multiple decks (R&D day vs quarterly), prefer the most recent and
note any change in assumptions across versions (a quiet change is a finding).

**When a document can't be fetched**: if a site blocks automated access or a paper is paywalled, don't try to get around it. Record the document as `blocked` and ask the user to supply the PDF if they have access.

**The text layer is triage, not extraction.** After locating the key pages by text, read them visually with `Read` (the visual pass in SKILL.md). Figures carry placebo-corrected CIs, error bars, and curve shapes the text layer omits; treat any number visible only in a figure as recoverable data, and cite it to the slide.

**When only the result is public** (paywalled paper, no deck): back out event
counts from the reported CI (SE of log HR = ln(upper/lower)/3.92, events ~
4/SE^2 at 1:1) and proceed to implied-effect mode. Label everything inferred.
Note the balanced-events approximation understates true event counts when the
death split is lopsided; if arm-level counts surface later, recompute with
Var = 1/d1 + 1/d2.

## Ingestion path B2: SEC filings (US-listed sponsors)

For any US-listed sponsor (subject or precedent), check EDGAR alongside the
sponsor's own pages, keeping requests low-volume. Full-text search the trial name, acronym, or NCT ID: JSON API
`https://efts.sec.gov/LATEST/search-index?q=%22<phrase>%22&startdt=YYYY-MM-DD&enddt=YYYY-MM-DD`
(UI at sec.gov/edgar/search). What each filing type yields:

- **8-K exhibits**: topline press releases and furnished decks.
- **S-1 / F-1**: for pre-commercial biotechs, often the most detailed public
  design disclosure outside the protocol (liability pushes toward
  completeness); check the Business section and pipeline discussion.
- **10-K / 10-Q**: design changes, enrollment status, readout timing.
  Compare consecutive filings: shifts in lawyer-vetted wording are signal;
  `assumptions_changed_quietly` applies to filings as well as decks.

Log SEC documents in `deck_manifest.json` like decks, with
`role: "sec_filing"` and the accession number or filing URL as the locator.

## Ingestion path C: Published paper

Methods section, "Statistical Analysis" or "Sample Size" paragraph. NEJM/JAMA/
Lancet are reliable about disclosing the original powering assumption. Also
pull the **observed** effect, CI, and observed control arm performance from
Results; these feed the post-readout (academic perspective) checks.

For protocol papers (e.g. "Design and rationale of the XYZ trial"), the entire
paper is the SAP; use it as the primary source.

## Ingestion path D: Manual

User supplies parameters directly. Echo them back in the table and tag
`source_type: user`.

## Endpoint classification heuristics

| Endpoint text contains | Classify as |
|---|---|
| ORR, response rate, proportion, % achieving, responder analysis, remission rate | `binary` |
| change from baseline, mean, score, LS-mean difference, ANCOVA | `continuous` |
| OS, PFS, EFS, DFS, RFS, TTP, time to, survival, KM, hazard ratio | `tte` |
| annualized rate, exacerbation, relapse rate, attack rate, events per patient-year, seizure frequency | `count` |
| mRS, ordinal shift, win ratio, composite | `ordinal`/`composite` (flag; use closest proxy) |

If the primary endpoint is a co-primary or hierarchical family, capture the
full hierarchy and the alpha split; the multiplicity handling is a Phase 5
red-flag check.

## Missing-parameter policy

The two parameters that are almost always findable: **n_total** (CT.gov) and
the **endpoint type** (CT.gov outcome text). With those alone you can compute
the implied effect size in Phase 2. Everything else is an audit input, not a
blocker. If `alpha` or `stated_power` is missing, use the conventional default
(two-sided 0.05, 80% for Phase 2, 90% for registrational Phase 3) and flag it
as an assumption in the memo.
