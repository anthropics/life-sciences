# Endpoint archetypes and indication conventions

The same statistical machinery powers an MCCB trial in schizophrenia and a
Schirmer test in dry eye, but the regulatory convention, comparator
behavior, and the trap that kills the trial differ. This file gives a
procedure that generalizes to any indication and a catalog that routes the
common cases fast.

## The procedure (this is what generalizes, not the catalog)

Before running any math in Phase 2, answer three questions and record the
answers in `design_params.json` under an `archetype` block:

1. **What is the regulatory endpoint convention in this indication?**
   Retrieve the FDA disease-specific guidance (search
   `site:fda.gov guidance [indication] developing drugs`) and the design of
   the 2-3 most recent approvals in the indication (labels, reviews, or
   registration trial papers). The convention is a design input, not a
   Phase 5 afterthought: it determines whether the design is single-primary,
   co-primary, sign-plus-symptom, replication-across-trials, or
   event-driven, and the powering math follows from that structure.
2. **What does the comparator arm do in this indication?** Placebo, vehicle,
   sham, active comparator, and natural history behave differently, and the
   comparator response is almost always the assumption sponsors get wrong.
   Pull it from the precedent trials, never from the sponsor's deck alone.
3. **Which archetype below is closest?** If none fits, the endpoint is
   probably simulation-sized (see "Other design shapes" in
   02-stats-methodology.md); find the design's published simulation
   assumptions and audit those.

If the indication has a qualified or de facto standard instrument (MCCB,
MADRS, ADAS-Cog, EASI, ACR core set, FEV1, mRS), the precedent search in
Phase 3 should be restricted to trials using that instrument; SDs and
placebo responses do not transfer across instruments.

## Archetype catalog

### A. Continuous symptom/function scale vs placebo

Psychiatry (MADRS, PANSS), pain (NRS), respiratory (FEV1), neurology
(motor scores). Mode: `continuous`. The killers are placebo response drift
(placebo arms improve more in recent trials than older ones; audit against
the last 5 years, not the classics) and SD inflation from Phase 2 to
Phase 3. MCID vs MDE matters for the payer perspective. Regional and site effects
are real in psychiatry; a US-only Phase 2 delta rarely survives a global
Phase 3.

### B. Co-primary across domains (cognition + function)

CIAS is the canonical case: FDA requires a cognitive primary (MCCB is the
de facto standard battery) AND a functional co-primary (interview-based or
performance-based functional capacity measure), because cognitive
improvement alone is not accepted as clinically meaningful. Alzheimer's
symptomatic trials carry the same structure (cognition plus global or
function). Mode: `continuous` per endpoint, then apply the co-primary
joint-power rule from 02-stats-methodology.md: per-endpoint 90% is 81% joint
at zero correlation, partially repaired by correlation. The audit
questions: is the stated power joint or per-endpoint, what correlation did
the sponsor assume between the co-primaries (it is usually undisclosed and
optimistic), and is the functional measure one FDA has actually accepted.
MCCB-specific traps: practice/learning effects shrink the placebo-adjusted
delta on repeat administration, and the functional co-primary typically has
the weaker signal, so it is the binding constraint; power the design on the
weaker endpoint, not the headline one.

### C. Sign plus symptom, replicated (ophthalmology surface disease)

Dry eye is the canonical case: FDA has historically required efficacy on
BOTH a sign (objective: staining, Schirmer tear production) and a symptom
(patient-reported: ocular discomfort, OSDI items), each supported by two
adequate trials, against vehicle. Mode: usually `continuous` or `binary`
per endpoint. The killer is vehicle response: the vehicle is itself a
lubricant, so vehicle arms improve substantially on both signs and
symptoms, and the placebo-adjusted delta is a fraction of the within-arm
change the deck headlines. Audit moves: power on the VEHICLE-adjusted
delta from precedent (not the sponsor's open-label or within-arm change),
check which specific sign endpoint the sponsor picked (Schirmer responder
cuts, staining subscales, and CAE-challenge designs have very different
precedent rates), and treat a sponsor switching sign endpoints between
trials as a multiplicity red flag. The replication structure also means the
relevant probability is winning the REQUIRED SET of trials, not one: two
independent trials each at 80% power is 64% joint, before correlation
arguments. Allergic conjunctivitis (CAE), blepharitis, and meibomian gland
disease inherit most of this archetype.

### D. Responder endpoints dichotomized from a scale

ACR20, PASI75/90, EASI75, clinical remission scores, IGA 0/1. Mode:
`binary`. Two traps. First, the placebo responder rate is exquisitely
sensitive to the baseline-severity entry criterion and the exact responder
definition; precedent must match both. Second, dichotomizing a continuous
measure discards information (asymptotic efficiency ~64% at a median split),
so a responder primary needs materially more patients than the underlying
continuous scale would; sponsors who power the responder endpoint using the
continuous effect size have understated N. Regression to the mean
guarantees placebo "response" whenever entry requires a minimum severity:
patients enroll at a high-water mark.

### E. Event-driven TTE with composite conventions

Oncology (OS, PFS) is the home case and is covered in
02-stats-methodology.md, including NPH and crossover. Cardiovascular
outcomes trials add the composite convention (MACE: CV death + MI +
stroke): power rides on the component mix, the softest component drives
event counts but the hardest carries the clinical story, and a composite
HR powered on total events can hide a null on mortality. Heart failure has
been moving to recurrent-event and win ratio primaries; route those to the
simulation-audit path. Audit move for any composite: get the precedent
component split and check whether the assumed HR is plausible for the
components that will actually contribute most of the events.

### F. Recurrent-event counts

COPD/asthma exacerbations, migraine days, HAE attacks, seizure frequency.
Mode: `count` (negative binomial). The named trap is overdispersion (k
understated), covered in 02-stats-methodology.md. Two indication-level
additions: entry criteria requiring a minimum event rate (e.g. >= 2
exacerbations last year) guarantee regression to the mean, so placebo arms
improve and the precedent placebo rate must come from trials with the SAME
entry threshold; and seasonal endpoints (asthma, migraine) make the
exposure window and enrollment calendar part of the rate assumption.

### G. Ordinal and global scales

Stroke (mRS shift), COVID severity scales, some derm IGA analyses. No
closed form in the script; follow the simulation-audit path in
02-stats-methodology.md and sanity-check with the nearest dichotomization
(mRS 0-2 vs 3-6 as `binary`), labeled as a proxy.

### H. Within-patient and paired designs

Ophthalmology fellow-eye, bilateral dermatology, crossover in stable
chronic disease. The unit of analysis is the pair; the effective variance
carries the within-patient correlation (crossover: 2*sd^2*(1-rho) per
period pair; see 02-stats-methodology.md for how to pass this to the script). The recurring error in published trials is analyzing eyes or
lesions as if independent, which overstates the effective N; if a sponsor's
power calc counts 2N eyes from N patients without a correlation adjustment,
flag it. Paired designs also cannot be rescued post hoc: if the correlation
assumption was wrong, the trial is simply mis-sized.

### I. Small-N rare disease

The z-approximations in power_calc.py are optimistic below ~30 per arm
(documented in the conventions section); rare disease trials often run
n=20-60 total. Expect exact tests, Bayesian designs, natural-history or
external controls, and within-patient designs. The audit shifts from "is
the power right" to "is the comparison interpretable": an external control
makes the effective alpha unknowable, and the right output is a labeled
approximation plus a statement of what the design cannot demonstrate.

### J. Prevention and vaccine efficacy

Case-driven, not patient-driven: VE = 1 - relative risk, and the
information unit is accrued cases, exactly analogous to events in TTE. The
control-arm assumption is the attack rate, which is epidemiology (season,
geography, variant) rather than trial precedent, and overestimating it is
the classic failure (trials extend because cases never accrue). Mode:
`binary` on case counts for a fixed-time read, or treat cases as events for
event-driven framing. Wittes 2002's noncompliance machinery applies with
full force here because prevention trials run for years.

## What this changes in each phase

Phase 1: record `archetype`, `regulatory_convention` (with the guidance or
approval citations), and `comparator_behavior` in `design_params.json`.
Phase 2: the archetype picks the mode and the structural correction
(joint power, replication probability, paired variance, dilution).
Phase 3: the archetype defines WHICH precedent parameters are load-bearing
(vehicle response for C, placebo responder rate at the same entry severity
for D, attack rate for J) and restricts the reference class to trials using
the same instrument and entry criteria.
Phase 5: archetype-specific red flags (sign-endpoint switching, co-primary
power stated per-endpoint, eyes counted as patients) go in the memo
alongside the generic ones.
