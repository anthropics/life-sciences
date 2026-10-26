# DSMB / interim monitoring perspective

## Core question

Given where the trial is right now, should it continue? This is the
mid-trial perspective: conditional power, futility, and sample size re-estimation.
It applies when the user has interim information (a disclosed interim look,
a blinded event count, an enrollment update) and wants to know what it
implies, or is designing the monitoring plan itself.

Unlike the other perspectives, the information state matters as much as the
design: what is known, at what information fraction, and under whose custody
(blinded vs unblinded).

## What to compute

1. **Conditional power.** P(final analysis significant | interim data), under
   three standard assumptions: the originally hypothesized effect, the
   current trend (interim estimate carried forward), and the null. Run
   `power_calc.py cp --z-interim <z> --info-frac <t> ...` for all three.
   Conventional read: CP under current trend below ~10-20% supports a
   futility stop; above ~80% the trial is on track.
2. **Efficacy boundary convention.** When the SAP names O'Brien-Fleming but
   not the exact spending function, different implementations give different
   interim critical values. State which
   convention you used and present the implied boundary as a range when the
   function is not disclosed, e.g. "boundary HR ~0.74-0.75 depending on the
   spending implementation." Reference values for a single interim plus
   final, two-sided 0.05: classic OBF gives z ~2.80 at 50% information, ~2.57
   at 60% and ~2.33 at 75%; Lan-DeMets OBF spending gives ~2.96, ~2.67 and
   ~2.34 in its one-sided 0.025 form, or ~2.77, ~2.53 and ~2.26 in its two-sided
   0.05 form.
3. **Futility boundary placement.** If the SAP discloses a futility rule
   (conditional power threshold, beta-spending, or a fixed HR boundary),
   check where the current estimate sits relative to it. If no rule is
   disclosed, note that nonbinding futility is the norm and DSMBs deviate.
4. **Information fraction tracking.** For TTE: events observed / events
   planned, not patients enrolled. Slow event accrual (control arm
   outperforming, or both arms doing well) pushes the readout right and is
   itself informative; compare blinded aggregate event rates against the
   design assumption.
5. **Sample size re-estimation read.** Blinded SSR on nuisance parameters
   (pooled event rate, pooled SD) is low-risk and common; an announced SSR
   tells you which assumption the sponsor got wrong, and the direction.
   Unblinded promising-zone SSR is a signal the interim effect landed in the
   marginal zone; treat an announced enrollment increase after an unblinded
   look as informative about effect size (likely smaller than hoped, not
   dead).
6. **Reading public DSMB actions.** Public DSMB actions carry information: a "continue without
   modification" after a planned efficacy look means the efficacy boundary
   was not crossed (bounds the interim effect from above); a trial that
   passes a futility look has CP above the futility threshold (bounds it
   from below). Quantify the implied effect range when a DSMB action is the
   only public information.

## Custody discipline

This perspective often runs on partial information. Be explicit about what is known
vs inferred, and never present an implied-effect range from a DSMB action as
an estimate; it is a bound. If the user is an actual trial team member with
unblinded data, the calculations are the same but flag that sharing
unblinded interim results outside the DSMB violates trial integrity.

## Output emphasis

A continue/stop dashboard: information fraction, conditional power under the
three assumptions, position vs any disclosed boundaries, and what each
plausible DSMB action would imply. When only public DSMB actions are known,
lead with the implied effect bounds.
