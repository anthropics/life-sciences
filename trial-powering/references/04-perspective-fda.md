# FDA reviewer perspective

## Core question

Is this an "adequate and well-controlled" investigation under 21 CFR 314.126,
and is the statistical plan consistent with ICH E9 / E9(R1)?

The agency does not require any specific power; it requires control of Type I
error and a pre-specified, interpretable analysis. Power is the sponsor's
problem. So this perspective focuses on alpha, multiplicity, estimand, and analysis
integrity rather than on whether the effect size assumption was rosy.

## What to check

1. **Type I error control.**
   - Overall alpha ≤ 0.05 two-sided (or 0.025 one-sided) for the primary
     family.
   - Multiplicity: if there are co-primaries, multiple doses, multiple
     populations, or a secondary hierarchy intended for label claims, is there
     a stated adjustment (Bonferroni, Hochberg, hierarchical/fixed-sequence,
     graphical)? Unadjusted multiple primaries is a red flag.
   - Interim analyses: alpha-spending function named, stopping boundaries
     pre-specified, DMC charter referenced.
2. **Estimand (ICH E9(R1)).** Is the target of inference defined? Population,
   variable, intercurrent-event handling (treatment policy / hypothetical /
   composite / while-on-treatment / principal stratum), and summary measure.
   Missing-data handling tied to the estimand (MI, MMRM, tipping point).
   Absence of an estimand statement in a post-2021 SAP is notable.
3. **Analysis population.** ITT (or FAS) for the primary superiority test;
   per-protocol as sensitivity. For NI, both ITT and PP should support the
   conclusion.
4. **Non-inferiority margin justification.** If NI: margin derived from
   historical evidence of the active control's effect (usually a fraction of
   the lower CI bound from placebo-controlled trials), constancy assumption
   discussed, assay sensitivity addressed.
5. **Endpoint acceptability.** Is the primary endpoint one FDA has accepted in
   this indication (check disease-specific guidance and recent approvals)? If
   it's a surrogate, is there an accelerated-approval precedent?
6. **Adaptive design compliance.** If adaptive: pre-specified adaptation
   rules, simulations supporting Type I error control, firewall between DMC
   and sponsor, consistent with the FDA Adaptive Designs guidance (2019).
7. **Single-study sufficiency.** If this is the only pivotal: does it meet the
   "single adequate and well-controlled study + confirmatory evidence" bar
   (large, multicenter, internally consistent, p well below 0.05, clinically
   meaningful, often mortality/irreversible-morbidity endpoint)?

## Relevant guidance to cite

Pull these via `WebSearch site:fda.gov` and cite the specific section:

- ICH E9 Statistical Principles for Clinical Trials
- ICH E9(R1) Addendum on Estimands and Sensitivity Analysis
- FDA Guidance: Multiple Endpoints in Clinical Trials (2022)
- FDA Guidance: Adaptive Designs for Clinical Trials (2019)
- ICH E20 Adaptive Designs for Clinical Trials (draft reached Step 2b June
  2025, consultation closed December 2025, finalization expected 2026; cite
  for confirmatory adaptive designs and check current status when used)
- FDA Guidance: Non-Inferiority Clinical Trials (2016)
- FDA Guidance: Demonstrating Substantial Evidence of Effectiveness (2019 draft / updates)
- Disease-specific guidance for the indication

## Output emphasis

A compliance-style checklist: each item, status (addressed / not addressed /
unclear), and a one-line note with the SAP/protocol citation. Then a short
paragraph on whether the design would likely survive a statistical review and
what the reviewer would push back on.
