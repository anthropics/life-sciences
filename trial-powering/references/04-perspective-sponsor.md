# Sponsor / operator perspective

## Core question

How should I size this trial, and how do I de-risk the design against the
assumptions I'm least sure about?

This perspective runs the same math forward instead of adversarially. The user is
choosing N, not auditing someone else's N.

## What to surface

1. **Sample size table.** N required for 80% and 90% power across a small grid
   of plausible effect sizes (e.g. optimistic / base / conservative), at the
   intended alpha and allocation. Include the dropout-inflated N. State how
   the target difference was chosen using the DELTA2 categories (anchor,
   distribution, opinion, review of evidence; see 02-stats-methodology.md);
   a Phase 2 point estimate without shrinkage is the weakest justification
   and reviewers know it.
2. **Cost of conservatism.** The N delta between powering for the Phase 2
   effect vs powering for a shrunk effect. Frame it as patients and (if the
   user supplies a per-patient cost) dollars and months.
3. **Design levers that buy power without adding N:**
   - Enrichment (biomarker-select, severity threshold) and its effect on the
     expected delta and variance
   - Covariate adjustment / ANCOVA (effective variance reduction via baseline
     `r^2`)
   - Allocation ratio (1:1 is most efficient for superiority; unequal makes
     sense for ethical/recruitment reasons but costs power per total N)
   - Endpoint choice (a more sensitive instrument, a longer follow-up, a
     responder cut that maximizes separation)
4. **Adaptive options.** Blinded sample size re-estimation (SSR) on nuisance
   parameters (SD, event rate, control response) is low regulatory risk and
   directly hedges the assumptions you're least sure of. Unblinded SSR /
   promising-zone designs are higher complexity; mention them if the
   uncertainty is on the effect size itself. Group-sequential with an interim
   for futility caps the downside.
5. **Event-driven vs calendar-driven (TTE).** If the primary is TTE, recommend
   event-driven final analysis so power is preserved even if accrual or event
   rate misses plan.

## Output emphasis

A decision table the operator can take to a steering committee: rows = effect
size scenarios, columns = N for 80% / N for 90% / months to readout / notes.
Follow with a short list of the two or three design changes that most cheaply
buy robustness.
