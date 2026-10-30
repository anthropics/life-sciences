# Phase 4 — Regulatory & trials research

**When:** `waypoints/biology_soc.json` exists, `waypoints/regulatory_trials.json` doesn't.

## Objective

Understand how to run clinical trials for this indication: what endpoints
regulators accept, what trials have shaped the field, and what failed.
This covers Sections 4 (Regulatory Path) and 5 (Key Trials).

## Actions

**Step 1: FDA/EMA accepted endpoints**
- Web search site:fda.gov for guidance documents specific to this indication
  (e.g., `"FDA guidance" "[indication]" clinical trial endpoints`)
- Web search site:ema.europa.eu for EMA guidelines and assessment reports for this indication
- Identify primary endpoints used in successful registrational trials
- Distinguish: clinical endpoints vs. surrogate endpoints vs. PROs
- Note any endpoints that are "reasonably likely to predict clinical benefit"
  (accelerated approval pathway)
- For indications without disease-specific guidance: note this, and check whether approval precedents exist before calling the path novel

**Step 2: Regulatory precedents**
- Web search site:fda.gov for approval packages related to this indication
- What endpoints led to approval? What trial designs were accepted?
- Any breakthrough therapy designations, fast track, or priority review precedents?
- Relevant advisory committee discussions or Complete Response Letters

**Step 3: Trial design parameters**
- Search the registry for Phase 3 trials in this indication to identify patterns: the connector's search or the
  public registry
- List every Phase 3 trial found (and every pivotal Phase 2 where there is no Phase 3), with
  its status. Open the registry record of each and take the primary endpoint from it. If there
  are more than about 15, do this for those started in the last ten years and those named in
  guidelines or among the landmark trials, and count the rest in section 6 under "Not followed up"
- For each one that is completed, or active with enrolment closed, run one web search for its
  results (trial name or NCT ID plus "results"). A registry record can lag a conference or
  press-release readout by months
- Typical trial sizes (enrollment targets)
- Typical trial durations (primary endpoint assessment timepoints)
- Common comparator arms (placebo, active comparator, SOC)
- Per-patient costs: run one literature search for published trial costs in this indication. Report only a figure that a cited source states; if none is found, say what was searched rather than estimating

**Step 4: Landmark trials**
- Identify the trials that shaped current SOC (usually 3-7; fewer if the field has fewer)
- For each: NCT ID, drug, sponsor, phase, key results, impact on practice
- Look each up in the registry and in the literature
- Take each landmark trial's result from its primary paper or its posted registry results, not
  from a review. If the literature tool does not return the primary paper, search the web for
  it; a subagent limited to one source family may step outside it for this. An identifier you
  remember (PMID, PMCID, DOI, NCT ID) may be fetched to check it, as URL_Validation allows, and
  may be cited once a tool has returned the record and its title and authors match
- If only a secondary source could be read, say so next to the result and in section 6
- For active sponsors in this indication, fetch their `/pipeline` or presentations
  page and read recent conference decks, following the PDF rules in the research standards
  (slide decks are figure-first)

**Step 5: Notable failures**
- Identify significant clinical trial failures in this indication
- What mechanism/approach was tested? Why did it fail?
- Lessons learned that inform future trial design
- Search the literature for review articles on failed approaches

**Step 6: Write `waypoints/regulatory_trials.json`**
Follow the schema from `references/waypoint-schemas.md`. Include `trial_landscape` counts.

## Guidance

- Regulatory path information varies enormously by indication maturity:
  - Well-established (e.g., MASH): disease-specific FDA guidance or clear approval precedents exist
  - Emerging (e.g., sarcopenia): limited or no formal guidance
  - Novel (e.g., ageing): no regulatory framework at all
- Be explicit about which category this indication falls into
- For landmark trials: prioritize trials that changed practice, not just
  the most recent ones. But a pivotal readout, approval, label change, withdrawal or
  Complete Response Letter from the last 24 months always belongs in the report (4.2, 5.1
  or 5.2), even if practice has not changed yet. Before writing the waypoint, run at least
  one web search for regulatory news and one for trial results in this indication, each
  limited to the last 24 months
- For failures: focus on mechanism-level lessons, not just "this drug didn't work"
- Three source families side by side, with subagents as SKILL.md sets out: web search site:fda.gov for guidance/approvals,
  the registry for trial patterns, the literature for reviews of trial history
