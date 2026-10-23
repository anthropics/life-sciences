# Payer / HTA perspective

## Core question

Even if the trial wins statistically, does the evidence package support
reimbursement at the price the sponsor will ask? HTA bodies (NICE, ICER,
IQWiG, G-BA, CADTH) and US payers re-run the powering logic with a different
target: not "p < 0.05" but "is the effect large enough, certain enough, and
in the right population to justify coverage."

This perspective also serves guideline committees and clinicians asking whether the
result is clinically meaningful; the analytic core (MCID vs powered effect)
is the same.

## What to check

1. **MCID vs powered effect.** Find the minimal clinically important
   difference for this endpoint and population (published MCID literature,
   prior HTA appraisals in the indication, FDA/EMA guidance language). If the
   trial is powered to detect an effect smaller than the MCID, it can win
   statistically and still fail appraisal. Report: MCID, powered effect,
   observed effect (post-readout), and which side of the MCID the CI's
   conservative bound sits on.
2. **CI lower bound, not the point estimate.** HTA cost-effectiveness models
   are stress-tested at the unfavorable bound. Compute the ICER-relevant
   quantity at the lower bound of the effect CI and note how much of the
   value case survives there. Germany's IQWiG applies an explicit version of
   this rule for continuous endpoints: it generally requires the CI of the
   standardized mean difference to clear 0.2 (not just zero) before crediting
   added benefit. Compute the expected SMD CI at the planned N and report the
   probability of clearing that threshold; trials powered at 90% against zero
   routinely have under 50% probability of clearing the IQWiG bar.
3. **Subgroup power.** Payers restrict coverage to subgroups where the
   evidence is strongest. Check which pre-specified subgroups have enough
   events/N for a stable estimate, and which restriction decisions the trial
   cannot statistically support either way. An underpowered subgroup cuts
   both ways: no evidence of benefit, no evidence of absence.
4. **Surrogate vs final endpoint.** If powered on a surrogate (PFS, response
   rate, biomarker), what is the surrogate-to-OS/QoL transfer evidence in
   this indication? Prior HTA appraisals often state explicitly how much
   credit they give the surrogate.
5. **Comparator relevance.** Powering against a comparator payers no longer
   consider standard of care weakens the package even if the trial wins.
   Check whether the control arm matches current treatment guidelines in the
   major markets.
6. **Indirect comparison feasibility.** If the relevant competitor was never
   trialed head-to-head, can a network meta-analysis or MAIC bridge it? That
   depends on shared comparators, endpoint definitions, and population
   overlap; note whether the trial design enables or blocks it.
7. **QoL / PRO powering.** Payers weight quality of life; PRO endpoints are
   usually secondaries with no formal powering. Note whether TTD/QoL
   analyses are alpha-protected or descriptive only.

## Output emphasis

A coverage-readiness read: effect vs MCID (with the CI bound), subgroups the
evidence can and cannot support, surrogate credit, comparator currency, and
the two or three things an appraisal committee will fight about. Cite prior
appraisals in the indication where they exist.
