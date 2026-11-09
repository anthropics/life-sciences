# Phase 2 — Research

## When

`current_phase == "research"`.

## Objective

Gather detailed information for each competitor, with subagents as SKILL.md
sets out. Focus depth on threat level — Level 1 gets
exhaustive coverage, Level 4 gets minimal.

## Target information

**Subject asset** (at Level 1 depth, for its row in the main table):

- Efficacy with SPECIFIC NUMBERS for each reported trial, with the
  population, line of therapy and comparator
- Safety with INCIDENCE RATES
- Dosing and route (dose, frequency, route of administration), as approved
  or as used in the latest trial
- Stage in the indication and region, and the next stated milestone

For each competitor, gather (depth varies by threat level):

**All competitors**:

- Mechanism of action (detailed)
- Development stage and timeline

**Landscape (once, not per competitor)**:

- What current SOC achieves, including the share of patients reaching
  treatment goals
- Unmet needs
- Approved and in-development combinations
- Why withdrawn or failed programs stopped

**Level 1–2 (marketed and late-stage)**:

- Efficacy with SPECIFIC NUMBERS (response rates, CI, p-values), with the
  population, line of therapy and comparator, for each pivotal trial
- Safety with INCIDENCE RATES (not just "common AEs")
- Dosing and route (dose, frequency, route of administration)
- Limitations and unmet needs
- Commercial context from public filings (if marketed)
- Head-to-head data (if available)

**Level 1 same-target**:

- How mechanism differs from subject asset
- Competitive threat assessment
- Timeline risk (who's ahead)

**Level 3–4 (early-stage)**:

- Mechanism differentiation
- Key differentiator if notable

## Stopping condition

- Subject asset: efficacy, safety and dosing and route researched, or those
  data confirmed not public
- Level 1 competitors: 100% researched with efficacy/safety data, or those
  data confirmed not public
- Level 2 competitors: ≥90% researched
- Level 3–4: At least identified with mechanism
- Or the budget in SKILL.md is reached. Then stop, record each competitor still
  short of its level under `gaps`, and go on to Phase 3 with
  `current_phase="synthesis"`, gaps and all. A map with stated gaps is worth
  more than research with no map
- No placeholder text (TBD, Unknown) for Level 1–2 efficacy/safety. Write
  "Not publicly available" only after checking the canonical source AND
  running at least one search aimed at disproving it, as the anti-fabrication
  protocol and Absence_Claims in `00-research-standards.md` require. For trial
  results, ClinicalTrials.gov often has none posted when a paper or abstract
  exists, so search PubMed, conference abstracts and sponsor press releases
  too. Only then does it count as researched, not as a gap

## Approach

If a subagent tool is available, use it for parallel research; running the same
research yourself is the fallback for when none exists. How you organize the work is up to you — by
category, by threat level, by data source, or however makes sense for this
landscape.

When spawning subagents:

- Give each one the content of `references/00-research-standards.md`
- Pass the ACTUAL list of competitors to research (not placeholders) and the
  subject asset profile, and give one subagent the subject asset itself
- Have them end each file with an `## Anomalies` section holding their
  `anomaly_flags`, and list any PDFs worth downloading
- Tell them which sources to use (see below)
- Have them write findings to markdown files in `<workdir>/waypoints/`
  (e.g. `research-level1-[target].md`). Without subagents, write the same files
  yourself
- Spawn multiple at once for parallel execution where the tooling allows

## Sources

- ChEMBL: Drug details, mechanism, target information
- Trial registry: Trial design, endpoints, results
- FDA documents (web search on accessdata.fda.gov, then fetch): Approved
  drug labels, safety information
- PubMed: Published efficacy and safety data
- Web search: Gap filling, conference data, recent news
- Licensed data, where the user has it: Pipeline and competitor data

## Outputs

1. Complete the research, including any subagents
2. Read the research files and consolidate into
   `waypoints/findings.json` (schema in `waypoint-schemas.md`)
   and carry every quote, source and anomaly across. Handle any PDFs that
   subagents listed, following Download_Safety
3. Check coverage — if gaps in Level 1–2, run targeted searches
4. **Checkpoint with user** before proceeding to Phase 3
5. Once the user continues, update `waypoints/progress.json`:
   `current_phase="synthesis"`. If gaps remain and the budget allows more
   research, stay in `"research"`; once the budget is reached, move to
   `"synthesis"` with the gaps recorded
