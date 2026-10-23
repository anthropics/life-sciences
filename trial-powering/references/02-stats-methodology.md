# Phase 2: Statistical methodology

The arithmetic lives in `assets/power_calc.py`; call it rather than doing
math in your head. This file ensures you call it with the right inputs and
interpret the output correctly.

## Core relationship

For a superiority test, power, alpha, sample size, and effect size are four
corners of one equation. Fix any three and the fourth is determined. Sponsors
publish N and (usually) alpha and power, so the **implied effect size** falls
out. That implied effect is the sponsor's bet, whether they say it out loud or
not.

Rule of thumb worth internalizing: required sample size scales roughly with
the inverse square of the standardized effect. Halve the true effect and you
need ~4x the patients; equivalently, if the true effect is half what the
sponsor assumed, power at the planned N collapses from ~90% toward ~30-40%.
This nonlinearity is why small assumption misses cause big readout misses.

## By endpoint type

### Binary (two proportions)

Inputs: `p_control`, `p_treat` (or `delta_abs`/`odds_ratio`), `n_per_arm`,
`alpha`, `sided`, allocation ratio.

`power_calc.py` uses the normal approximation with unpooled variance for power
and pooled variance for the test statistic (matches most SAPs). For ORR-style
single-arm vs historical control, use `--single-arm --p_null <historical>`.

Watch for: continuity correction (rarely used in modern SAPs), exact vs
asymptotic (matters when N per arm < ~30 or p near 0/1).

### Continuous (two means)

Inputs: `delta` (mean difference), `sd` (common or per-arm), `n_per_arm`,
`alpha`, `sided`, optional `ancova_r2`.

The standardized effect (Cohen's d = delta / sd) is what drives power. If the
sponsor reports delta but not SD, the implied-effect calculation returns d
directly. If an ANCOVA with baseline covariate is planned, effective variance
is `sd^2 * (1 - r2)`; sponsors sometimes bake this in without saying so, which
makes their N look small.

Watch for: SD taken from a Phase 2 that was smaller, shorter, or in a
healthier population (SD usually inflates in Phase 3).

Dropout inflation cuts both ways on longitudinal endpoints. The standard
n/(1-q) inflation is conservative when the analysis is MMRM, because
dropouts contribute partial information from their earlier visits; sponsors
sizing this way have a hidden cushion. But under a treatment-policy estimand
with retrieved-dropout imputation (now standard in obesity and diabetes),
the estimated effect itself is diluted by off-treatment data, so the honest
adjustment is to shrink the effect (like the noncompliance section), not
just inflate N. A sponsor who inflates N for dropout while powering on the
on-treatment effect under a treatment-policy estimand has mismatched the
power calc to the analysis; flag it.

### Time-to-event (log-rank / Cox)

The information unit is **events**, not patients. Schoenfeld's formula:

```
events_required = (z_{1-α/s} + z_{1-β})^2 / (log HR)^2 * (r+1)^2 / r
```

where `r` is the allocation ratio (treat:control). `power_calc.py
--type tte` computes events↔HR↔power and, given accrual + follow-up +
control median, the N needed to observe those events.

The implied HR from a stated event target is the headline number. If the
sponsor only states N (not events), back out the expected event count from
their control median, accrual, and follow-up assumptions, then compute the
implied HR; flag that the event assumption is itself a sensitivity axis.

Watch for: non-proportional hazards (delayed separation in IO, crossing
curves), informative censoring, interim looks consuming alpha (OBF spends
little early; Pocock spends a lot).

**Non-proportional hazards, quantified.** The log-rank statistic tracks the
event-weighted average log HR. Under delayed separation (typical of
checkpoint inhibitors), a fraction f of events landing before separation
carries no effect, shrinking the effective log HR by ~(1-f) and inflating
required events by ~1/(1-f)^2. Worked example: separation at 3 months with
control median 6 months puts ~30% of events pre-separation, turning a true
post-separation HR of 0.70 into an average near 0.78 and roughly doubling
the events needed. Audit moves: estimate f from precedent KM curves,
recompute power at the diluted average HR, and check whether the SAP names
an NPH-robust analysis (weighted log-rank, MaxCombo, or RMST). Powering a
delayed-separation mechanism on plain log-rank at the post-separation HR is
the most common quantitative powering error in modern IO trials. Crossover
into the control arm at progression is the OS mirror image: it is drop-in
(noncompliance section) and attenuates the OS HR toward 1 by roughly the
crossover fraction times the effect, which is why OS often misses where PFS hit.
Check crossover provisions and precedent crossover rates before believing
an OS power claim.

### Co-primary endpoints

Winning on ALL of K endpoints (co-primaries) means each is tested at full
alpha, but Type II error compounds: joint power is the product of the
marginals at zero correlation, and positive correlation helps only partly.
Two co-primaries each at 90% give 81% joint at zero correlation and about 83%
at a correlation of 0.5 (Sozu et al.; FDA Multiple Endpoints guidance 2022). The audit question: is the stated power
JOINT or PER-ENDPOINT? "90% power" on two co-primaries almost always means
per-endpoint (~81-85% joint at realistic correlation). Winning on ANY of K
is the opposite trade: alpha must be split or gatekept, and each endpoint's
powered effect recomputed at its allocated alpha. `power_calc.py` handles
endpoints marginally; multiply marginal powers for a conservative joint
estimate (a lower bound when the endpoints are positively correlated) and say
so in the memo.

### Count / rate (negative binomial)

Inputs: `rate_control`, `rate_ratio` (or `rate_treat`), exposure time,
overdispersion `k`, `n_per_arm`.

Variance = `rate * (1 + rate * k * t) / t`. Overdispersion is the parameter
sponsors most often understate; MS, COPD, HAE, and migraine trials routinely
show k well above what early-phase data suggested. `power_calc.py --type
count` implements the Keene/NB approximation.

### Non-inferiority

Same machinery, but the null is `delta = -margin` (or `HR = margin`) instead
of zero/one. Always one-sided; convention is alpha=0.025. The choice of margin
is a regulatory question (Phase 5 / FDA perspective), not a power question, but the
power calc must use the margin, not zero.

Two NI-specific sizing traps. First, analysis population: in NI trials ITT is
anti-conservative (sloppy conduct pulls both arms together and fakes
equivalence), so regulators typically require NI to hold in both ITT and
per-protocol (or micro-ITT for anti-infectives); size to the evaluable
population and gross up enrollment by the precedent evaluable-to-randomized
yield, which runs well below 1 in anti-infective trials. Second, a lower true
event/cure rate moves variance, usually upward (p toward 0.5 inflates pq), so
"rates might be worse than assumed" generally means MORE patients, and it
also stresses constancy: the margin was justified at the historical rates,
not the drifted ones.

## Choosing the target difference (DELTA2)

The target difference is a choice, not a measurement. DELTA2 (Cook et al.,
BMJ 2018) names the legitimate anchors: anchor methods (patient-meaningful
change on another measure), distribution methods (fractions of SD), opinion
seeking, and review of evidence (prior trials, meta-analysis). When
auditing, ask which of these the sponsor's delta traces to; "our Phase 2
point estimate" is review-of-evidence with no shrinkage, the weakest form.
When sizing as a sponsor, document the method; reviewers expect it.

## Assurance (probability of success) vs power

Power conditions on one assumed effect. Assurance averages power over a
prior on the effect (O'Hagan; Bayesian expected power) and is the standard
industry PoS metric. The Phase 3 audit produces exactly its inputs: run
`power_calc.py assurance --prior-mean <effect> --prior-sd <spread>` with the
precedent-informed center and reference-class spread, and report assurance
next to point power. The gap is itself a finding: 90% power but 60%
assurance means success depends on the point assumption being right. Two
cautions: assurance is capped well below 1 even for huge trials (prior mass
near the null never converts), and the prior must come from the audited
reference class, not the sponsor's deck.

## Noncompliance, drop-in, and dilution (intention-to-treat sizing)

Noncompliance dilutes the effect rather than removing patients. Under ITT,
dilute BOTH arms when sizing: treated rate = (1-t)*p_treat + t*p_control
(drop-out fraction t), control rate = (1-c)*p_control + c*p_treat (drop-in
fraction c). Required N inflates by ~1/(1-c-t)^2 (Wittes, Epidemiol Rev
2002;24:39-53): 20% drop-out + 5% drop-in costs ~1.78x the sample.

Sizing on the per-protocol subset instead is a trap worth naming whenever a
sponsor or operator proposes it: excluding non-compliers preserves nominal
power on paper but breaks randomization, and the resulting estimate is
biased. High per-protocol power is deceptive, not reassuring.

## Approximation boundary: "X-year event rate" designs

"Event within K years" plus staggered recruitment is really time-to-event
with variable follow-up. Two-proportion formulas understate N vs
survival-based methods (Lachin-Foulkes, Lakatos) because late enrollees
contribute less information; in Wittes' worked example the gap is ~30-40%
(3,300 vs ~2,400). Use tte mode with
`--n-from-events --med-ctrl --accrual --followup` (a simplified
Lachin-Foulkes) whenever the pattern appears, and present the binary number
only as a labeled lower bound.

## Other design shapes and their handling

- **Competing risks.** When a competing event (death in a non-mortality
  endpoint) censors patients, person-time is lost and the effective event
  count drops below plan. The fix at audit time: discount expected events by
  the competing-event rate from precedent and rerun the grid. If the
  competing risk plausibly differs by arm, flag it; no sample size
  adjustment fixes differential informative censoring.
- **Ordinal, win ratio, hierarchical composites** (mRS shift, Finkelstein-
  Schoenfeld, win odds). No closed form here. Audit approach: find the
  design's published simulation assumptions (these designs are almost always
  simulation-sized and the assumptions are usually in the SAP or design
  paper), sanity-check the implied per-component effects against precedent,
  and run the nearest analytic proxy (binary on the dominant component, or
  tte) clearly labeled as a proxy.
- **Cluster randomization.** Inflate N by the design effect
  1 + (m-1)*ICC for cluster size m; trials that ignore this are underpowered
  by construction. Rare in drug registration, common in health-services and
  vaccine field trials.
- **Crossover designs.** Within-patient comparison; effective variance is
  2*sd^2*(1-rho) per period pair. The continuous machinery is a two-sample
  formula, so for N pairs pass `--n1 N` and `--sd <SD of differences>/sqrt(2)`.
  If the paired SD isn't disclosed, flag that the powering is opaque.

## Numerical conventions in power_calc.py

Stated here so memos can answer "which test does this assume": binary uses
pooled variance under H0 (the standard SAP z-test; unpooled-H0 tools differ
by ~1pp at unequal allocation), continuous uses the normal approximation
(optimistic by ~1 patient/arm below n=30/arm; use exact-t software for small
trials), the NI branch evaluates variances at the true rates (Blackwelder),
and `grid --vary-dropout` models attrition (information loss), not ITT
dilution, which belongs in the effect parameters per the noncompliance
section above.

## Interim analyses and alpha spending

If the design has interims, the final-analysis critical value is larger than
the nominal z, so power at the final look is slightly lower than the naive
calc suggests. `power_calc.py` accepts `--alpha-spent <x>` to use the
remaining alpha at final, which is conservative; if the SAP gives the final
nominal alpha, pass that as `--alpha` instead. If the spending function is OBF and there is one
interim at ~50% info, the effect is small (final alpha ~0.048); for Pocock or
many looks it matters more. If the SAP isn't disclosed, run both a naive and a
"one OBF interim" version and report the range.

Group-sequential designs can also **stop early for efficacy**, which raises
overall power above the per-look power. Don't try to model this precisely; if
relevant, note it qualitatively.

## What to write in `power_recon.json`

```json
{
  "forward_check": {
    "inputs_used": {...},
    "computed_power": 0.89,
    "sponsor_stated_power": 0.90,
    "reconciles": true,
    "notes": "matches within rounding"
  },
  "implied_effect": {
    "metric": "HR" | "abs_risk_diff" | "cohens_d" | "rate_ratio",
    "value": 0.72,
    "at_alpha": 0.05, "at_power": 0.90, "at_n": 600,
    "interpretation": "Trial is sized to detect a 28% hazard reduction"
  },
  "mde_at_80": {"metric": "...", "value": ...},
  "break_even_effect": {"metric": "...", "value": ..., "note": "true effect at which power = 50%"}
}
```

If the forward check does not reconcile (computed power differs from stated by
more than ~3 absolute points), list the candidate explanations: undisclosed
dropout inflation, ANCOVA adjustment, one-sided vs two-sided mismatch, interim
spend, or a different test (e.g. stratified vs unstratified). This discrepancy
goes in the memo as a finding.
