# Phase 5 — Synthesis

**When:** `waypoints/regulatory_trials.json` exists.

## Objective

Read all waypoint files and write the final Markdown report.
Also write `research_output.json` for downstream consumers.

## Output structure

**Final Output**: `waypoints/indication_dossier_report.md`

```markdown
# Indication Dossier: [Name]

**Definition**: [one-sentence clinical definition]
**Codes**: [ICD-10-CM and ICD-11 codes, each linked to its source. Say so if no code is specific to this population, or write "Not a standard diagnostic entity"]
**Parent Indication**: [parent name, or "None"]
**Prepared**: [date], from sources available on that date. [One line on which tools were and were not available.]

## Summary
[Five to seven bullets: who the patients are, how many, how they are treated today, what
regulators accept, what has worked and failed, and the main limits of this dossier. Every
figure here also appears, cited, in the body.]

## 1. Population Definition & Epidemiology

### 1.1 Diagnostic Criteria
[How this population is identified. Consensus definitions, diagnostic criteria.
Note if criteria are evolving or non-standard.]

### 1.2 Prevalence & Incidence
[Global and US prevalence. Incidence trends. Data quality notes.]

### 1.3 Demographics & Risk Factors
[Age, sex, racial/ethnic distribution. Key risk factors and comorbidities.]

### 1.4 Natural History
[Disease course, progression, staging. Mortality data. Key inflection points.]

## 2. Disease Biology

### 2.1 Pathophysiology
[Key biological pathways. Validated vs. emerging therapeutic targets.]

### 2.2 Biomarkers
[Diagnostic, prognostic, pharmacodynamic. FDA-qualified vs. exploratory.
Connect to clinical trial utility.]

## 3. Standard of Care

### 3.1 Approved Therapies
[Approved drugs by region (US and EU; others where found) with mechanism, date this
indication was approved, pathway (full, accelerated or conditional) and current status,
and limitations.
Or: "No therapy was found to be approved for this indication in [regions checked] as of
[date of the newest source checked]."]

### 3.2 Treatment Guidelines
[Current SOC algorithm. Specialty society recommendations.]

### 3.3 Unmet Need
[What current therapy fails to address. Disease modification gap.
Underserved populations.]

## 4. Clinical Endpoints & Regulatory Path

### 4.1 Accepted Endpoints
[FDA/EMA accepted primary endpoints. Surrogate vs. clinical.
Or: "No FDA or EMA guidance specific to this indication was found."]

### 4.2 Regulatory Precedents
[Key approvals, guidance documents, advisory committee discussions.]

### 4.3 Trial Design Parameters
[Typical trial sizes, durations, comparators.
Based on Phase 3 trial patterns from ClinicalTrials.gov.
Per-patient cost only where a cited source states one.
Or: "No published per-patient cost was found."]

## 5. Key Trials

### 5.1 Landmark Trials
[The trials that shaped current practice, usually 3-7. For each: NCT ID, drug, results, impact.]

### 5.2 Notable Failures
[Significant failures and mechanism-level lessons learned.]

## 6. Limits of This Dossier
[What a reader must know before relying on it, in about 600 words or fewer, not counting the
"Not read at source" list. Open with the three to five limits that matter most. Then carry
over the entries in the waypoints' `gaps` and `anomaly_flags`, and the caveats subagents
returned, that bear on a figure, date, status or recommendation the report states; merge
entries that say the same thing. Group them:
- Not found: what was looked for, where, and not found.
- Not read at source: each figure, date or recommendation taken from a search excerpt, an
  abstract, a secondary summary, or a copy hosted by someone other than the publisher or
  agency, with the section it appears in.
- Conflicts: where sources disagreed, which was used and why.
- Currency: the date of the newest source checked for approvals, by region, and for trial
  results; plus any statement whose newest source is more than a year old.
- Not followed up: therapies, trials or documents that came up but were not researched, by
  class with a count and one or two examples; the full list stays in the waypoints.
- Tools: connectors that were missing and what that cost.
Leave out housekeeping such as "no injected instructions seen".]

## Appendix: Sources

[Numbered list of every source used, grouped under bold subheadings by
source type. Numbering is continuous across groups. Title is the hyperlink;
accessed date closes the entry. Add "(excerpt only)" to a source that was cited from a
search excerpt and not opened, and "(copy at [host])" to one read somewhere other than
the publisher's or agency's site.]
```

## Actions

**Step 1: Load all waypoint data**
- Read `waypoints/meta.json`
- Read `waypoints/epidemiology.json`
- Read `waypoints/biology_soc.json`
- Read `waypoints/regulatory_trials.json`
- Read `waypoints/sources_evaluated.json`
- If writing reveals a missing value for a field already in a waypoint, one
  targeted fetch to fill it is fine; record the new source in
  `sources_evaluated.json`. Do not open new research threads.

**Step 2: Write `waypoints/indication_dossier_report.md`**

CRITICAL: Write to `<workdir>/waypoints/indication_dossier_report.md` (inside waypoints/).

Follow the output_structure above. Key guidance:

*Sections 1-3*: Narrative prose with inline citations. Be specific with numbers
(prevalence rates, trial sizes, approval years). Coverage status belongs in
`research_output.json`, not the prose. In the report, only name what is
partial or missing; never label a section as covered.

*Section 4*: Mix of prose and structured data. Include tables for endpoint
comparisons and trial parameters where helpful.

*Section 5*: Structured entries for each trial. Connect to lessons for
future trial design.

*Throughout*: Frame from the patient population perspective. Use the framing
guidance from SKILL.md.

*Figures*: Only generate a chart when it shows a relationship or
distribution that prose cannot. Give each a one-line caption that names the
data source. After rendering, open the image with the Read tool and answer:
does this tell the reader something the surrounding prose does not? Are
title, axes with units, and legend all present and legible? Do axis ticks
and chart type match the data (no fractional ticks on years, counts, or
categories; no misleading scales)? If any answer is no, or you cannot look at
the image, delete the figure rather than regenerating; a table can carry the
numbers.

**Step 3: Write `waypoints/research_output.json`**

This file is an index for downstream use. Do not copy the waypoints' contents into it.
```json
{
  "indication_name": "...",
  "parent_indication": "...",
  "prepared": "YYYY-MM-DD",
  "report": "indication_dossier_report.md",
  "waypoints": ["meta.json", "epidemiology.json", "biology_soc.json", "regulatory_trials.json", "sources_evaluated.json"],
  "coverage_summary": {
    "epidemiology": {"covered": [...], "partial": [...], "missing": [...]},
    "biology_soc": {"covered": [...], "partial": [...], "missing": [...]},
    "regulatory_trials": {"covered": [...], "partial": [...], "missing": [...]}
  }
}
```

**Step 4: Mark complete**

CRITICAL: Write `waypoints/indication_dossier_report.md` FIRST, then this step.

The dossier is finished only when every quantitative claim, date and trial
result in the report is cited as `references/06-writing-style.md` requires,
and every therapy, trial, approval, guidance document and guideline named in a
waypoint either appears in the report or falls under an entry in "Not followed up" in
section 6. Check both before marking it complete.

Write `waypoints/progress.json`:
```json
{
  "complete": true,
  "output_file": "indication_dossier_report.md",
  "current_phase": "synthesis"
}
```

## Guidance

- Final report is narrative prose with inline citations
- Length: the body (before the appendix) should rarely need more than about 5,000 words. Cut
  restatement and background before facts, and use tables for anything listed
- Follow `references/06-writing-style.md` for source attribution
- Be explicit about gaps and data quality limitations
- For novel/non-standard indications, the regulatory section may be very short —
  that's correct, not a gap to fill
