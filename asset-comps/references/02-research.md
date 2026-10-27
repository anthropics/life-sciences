# Phase 2 — Research

## When

`current_phase == "research"`.

## Objective

Gather detailed information for each competitor, with subagents as SKILL.md
sets out. Focus depth on threat level — Level 1 gets
exhaustive coverage, Level 4 gets minimal.

## Target information

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

- Level 1 competitors: 100% researched with efficacy/safety data, or those
  data confirmed not public
- Level 2 competitors: ≥90% researched
- Level 3–4: At least identified with mechanism
- Or the budget in SKILL.md is reached. Then stop, record each competitor still
  short of its level under `gaps`, and go on to Phase 3. A map with stated gaps
  is worth more than research with no map
- No placeholder text (TBD, Unknown) for Level 1–2 efficacy/safety. Where
  the canonical source has been checked and nothing is disclosed, write "Not
  publicly available": that counts as researched, not as a gap

## Approach

If a subagent tool is available, use it for parallel research; running the same
research yourself is the fallback for when none exists. How you organize the work is up to you — by
category, by threat level, by data source, or however makes sense for this
landscape.

When spawning subagents:

- Give each one the content of `references/00-research-standards.md`
- Pass the ACTUAL list of competitors to research (not placeholders) and the
  subject asset profile
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
   `current_phase="synthesis"`. If gaps remain, stay in `"research"`
