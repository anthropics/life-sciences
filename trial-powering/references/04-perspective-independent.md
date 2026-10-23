# Independent reviewer perspective (and academic post-readout perspective)

## Core question

Is this trial set up to succeed on realistic assumptions, or does the sponsor
need everything to break their way? What has to be true for the primary
endpoint to hit, and how likely is that?

An independent reviewer reads a powering analysis critically. Company decks
present a sponsor's case. The job is to find where the assumptions are
fragile.

## Pre-readout: what to surface

1. **The implied delta, in plain language.** "Powered for HR 0.72; class
   median in this setting is HR 0.80-0.85." One sentence, with the precedent
   comparison inline.
2. **Which single assumption is doing the most work.** Usually it's either the
   effect size or the control arm. Name it.
3. **The "what has to be true" table.** For each assumption: sponsor value,
   precedent median, verdict, and the power if precedent is right instead.
4. **The break-even line.** "If the true HR is worse than 0.79, this is a coin
   flip."
   When the Phase 3 audit yields a precedent spread, also report assurance
   (`power_calc.py assurance`): power averaged over the precedent prior.
   "90% power but 64% assurance" is the fragility story in one number.
5. **Early signals.** Things that update probability without waiting for
   data: enrollment pace vs plan, blinded event rate disclosures, DSMB
   continuations, sample size re-estimation triggers, protocol amendments that
   touch the primary endpoint or stat plan.

## Post-readout (academic perspective shares this section)

When the result is already out, the questions flip:

- **Was the observed effect within the powered range?** If the trial "just
  missed", was it underpowered for the effect that actually showed up, or did
  the drug simply not work?
- **Control arm reality check.** Did the control arm land where the sponsor
  assumed? If it overperformed, quantify how much of the miss that explains.
- **Fragility.** For a positive binary/tte result: how many events flipping
  the other way would erase significance (fragility index)? For continuous:
  how sensitive is the p-value to the SD used?
- **Post-hoc power is not informative on its own.** Don't report "observed
  power" as if it means something; it's a transform of the p-value. Instead
  report the CI on the effect and whether the powered-for effect sits inside
  it.
- **Compute realized information from the CI, not the stated event count.**
  SE of log HR = ln(upper/lower)/3.92. If the CI-implied information
  disagrees materially with the reported events (lopsided event split,
  heavy stratification, early censoring), flag the gap; it changes every
  downstream power statement and is sometimes itself a finding.
- **Multiplicity / selective reporting.** Did the topline match the
  pre-specified hierarchy? Any endpoint or population swap between SAP and
  press release?

## Output emphasis

Lead with the verdict, then the table, then the sensitivity headline. Keep the
mechanics (formulas, alpha spending) in an appendix. The reader wants to know
whether the design holds up; give them the one-paragraph answer first.
