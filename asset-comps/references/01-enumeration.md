# Phase 1 — Enumeration

## When

`progress.json` doesn't exist OR `current_phase == "enumeration"`.

## Objective

Find ALL potential competitors for the subject asset. Be exhaustive — it's
better to include too many than to miss a key competitor.

## Subject asset

Start by profiling the subject asset from sources: company, target, mechanism,
stage, key clinical data and any stated development timeline. Tiering and
every later comparison depend on this, so don't fill it in from memory. If it
can't be identified, stop and ask the user.

## Target information

For each competitor, capture:

- Drug name (generic and brand if applicable)
- Company/sponsor
- Target and mechanism class
- Development stage, with its source
- Category (`same_target`, `target_specific_soc`, `baseline_soc`,
  `adjacent_mechanism`, `other_pipeline`)

## Stopping condition

**Search until saturation.** You're done when additional searches stop
yielding new competitors. Specifically:

1. You've tried multiple search strategies:
   - Trial registry search by indication (paged through all results)
   - ChEMBL target and drug search by subject's target AND related pathway
     targets
   - Web search for competitive landscape overviews
   - Web search for treatment guidelines (to find SOC)
   - Sponsor searches for known players in the space

2. Cross-validation searches return no new names:
   - "[indication] competitive landscape [current year]" → no new drugs
   - "[target] pipeline" → no new programs
   - Names from web results already in your list

3. You've captured competitors across the full landscape:
   - Same target (direct competitors)
   - Same mechanism class (target-specific SOC)
   - Different mechanisms, same patients (baseline SOC, adjacent)
   - Pipeline across stages (Phase I/II/III)
   - Withdrawn, discontinued or failed programs (`withdrawn_failed`)

If you find new competitors, keep searching. If searches return only drugs
you've already listed, you've reached saturation.

## Tools guidance

**Licensed data**, where the user has it, is the preferred starting point for
pipeline and competitor data.

**Trial registry** (primary source):

- Search by indication to get the indication-wide pipeline
- Filter by intervention for SOC drugs
- Search by sponsor for known players' programs in the space
- Page through results to get complete lists

**ChEMBL** (target/mechanism enumeration):

- Target search by subject's target AND pathway targets
- Drug and mechanism lookups to identify same-target programs and mechanism
  class
- Compound search for known SOC molecules

**Web search** (gap filling):

- "[indication] discontinued program" and "[target] trial terminated" for
  withdrawn and failed programs
- "[indication] competitive landscape" for overview
- "[indication] treatment guidelines" for SOC identification
- "[target] inhibitor pipeline" for same-target programs
- "[indication] phase 3 pipeline [current year]" for recent entrants

## Threat levels

Assign threat level to each competitor (1=highest priority for research). If a
competitor fits more than one level, use the highest one its target, mechanism
class or stage reaches, but "ahead of the subject asset" applies only to
programs with the same target or mechanism class:

**Level 1** (direct threats — must research exhaustively):

- Same target programs, at any stage
- Target-specific SOC (marketed competitors with same mechanism)
- Phase III programs in indication
- Same target or mechanism class programs ahead of the subject asset

**Level 2** (important context — moderate research):

- Phase II programs
- Adjacent mechanism programs (Phase II+)
- Baseline SOC therapies

**Level 3** (awareness — brief coverage):

- Phase I programs (non-same-target)
- Early adjacent mechanism programs

**Level 4** (catalog only — list but minimal research):

- Preclinical programs (other than same target)
- Withdrawn or failed programs, unless their failure bears directly on the
  subject asset (then Level 2)

## Outputs

1. Create `<workdir>/waypoints` directory
2. Write `waypoints/enumeration.json` with all competitors (schema in
   `waypoint-schemas.md`)
3. **Checkpoint with user** before proceeding to Phase 2
4. Once the user continues, update `waypoints/progress.json`:
   `current_phase="research"`
