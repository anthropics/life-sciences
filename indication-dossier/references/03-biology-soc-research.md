# Phase 3 — Biology & standard-of-care research

**When:** `waypoints/epidemiology.json` exists, `waypoints/biology_soc.json` doesn't.

## Objective

Understand the biology underlying this patient population and how they're
currently treated. This covers Sections 2 (Disease Biology) and 3 (Standard of Care).

## Actions

**Step 1: Pathophysiology**
- Search the literature for recent review articles on the disease mechanism
- Identify key biological pathways involved
- Note which pathways are validated drug targets vs. emerging hypotheses
- For complex multi-pathway conditions, focus on the pathways most relevant
  to therapeutic intervention

**Step 2: Biomarkers**
- Diagnostic biomarkers: how is the condition confirmed?
- Prognostic biomarkers: what predicts disease progression?
- Pharmacodynamic biomarkers: what can measure drug effect?
- Note each biomarker's regulatory standing: formally qualified (check the public lists of
  FDA's Biomarker Qualification and Drug Development Tool programs, and EMA's qualification
  opinions), accepted in guidance or by approval precedent, cleared or approved as a test,
  or exploratory
- This directly informs Section 4 (endpoints)

**Step 3: Approved therapies**
- Web search site:fda.gov for approved drugs for this indication
- Web search site:ema.europa.eu for medicines authorised in the EU for this indication; note
  approvals in other regions when they come up
- For each approved therapy: mechanism, approval year, key limitations
- Confirm each approval year against the canonical source named in
  `references/00-research-standards.md` (Drugs@FDA, or FDA's CBER product lists). A drug
  first approved for another condition may have gained this indication later, so
  record the year this indication was approved
- Note if therapies are approved for this specific indication or used off-label; the
  indications section of the current label settles it, and the label also shows any
  boxed warning
- Identify major drug classes and their positioning (first-line, second-line, etc.)
- For indications with no approved therapies: state this explicitly

**Step 4: Treatment guidelines**
- Web search for specialty society guidelines (NCCN, AASLD, ATS/ERS, AGS, etc.)
- If a guideline's own text cannot be opened, use the society's summary or a peer-reviewed
  article that quotes it, and say in section 6 that the guideline itself was not read. Do not
  cite a drug's promotional site for what a guideline recommends
- Current standard of care algorithm
- How guidelines differ between regions (US vs EU) if relevant
- Recent guideline changes and their implications

**Step 5: Unmet need**
- What does current therapy fail to address?
- Patient populations underserved by existing treatments
- Disease modification vs. symptom management gap
- Quality of life impact not addressed by current therapies

**Step 6: Write `waypoints/biology_soc.json`**
Follow the schema from `references/waypoint-schemas.md`.

## Guidance

- For pathophysiology: aim for a working understanding — key pathways and
  why they matter for therapy, not textbook depth
- For biomarkers: explicitly connect to clinical trial utility
- For approved therapies: focus on limitations — this is what defines the unmet need
- Three source families side by side, with subagents as SKILL.md sets out: literature for biology, web search for
  guidelines, web search site:fda.gov and site:ema.europa.eu for approvals
