# Phase 5: Red-flag checklist

Screen the design against this list. For each flag that fires, work the
finding into the memo in plain English with its evidence ("the control-arm
assumption is anchored on seven-year-old data; contemporary trials run 2-3
months better"). The snake_case flag ids exist only for
`powering_output.json`; they never appear in the memo. Don't mention flags
that clearly don't apply.

## Powering assumption flags

| Flag id | Fires when | Why it matters |
|---|---|---|
| `effect_vs_class` | Implied effect is more favorable than the class median by a material margin | The single biggest predictor of Phase 3 misses |
| `phase2_anchoring` | Effect assumption equals the asset's open-label or single-arm Phase 2 result with no shrinkage | Phase 2 effects shrink in randomized Phase 3; powering at the unshrunk number is optimistic |
| `control_arm_stale` | Control arm assumption is anchored on data >5 years old or from a different line/population | Contemporary control arms outperform historical ones; this is the most common "drug worked but trial missed" story |
| `sd_from_small_n` | SD (or dispersion) taken from a Phase 2 with N<50 | Variance is underestimated in small samples; Phase 3 SD is usually higher |
| `dropout_low` | Assumed dropout below precedent median for trials of similar length/burden | Understated dropout inflates effective N on paper |
| `events_optimistic` (TTE) | Event accrual assumption requires a control median worse than precedent | Slow events delay readout and may force calendar-driven analysis at lower power |

## Design / analysis flags

| Flag id | Fires when | Why it matters |
|---|---|---|
| `no_multiplicity` | Multiple primaries / doses / populations with no stated alpha control | Inflated Type I error; FDA will not accept nominal p-values |
| `one_sided_unjustified` | One-sided alpha=0.05 in a superiority setting without rationale | Effectively doubles alpha vs convention; uncommon outside NI |
| `endpoint_changed` | Primary endpoint or analysis population changed after enrollment started (CT.gov history) | Possible outcome-driven amendment; credibility hit |
| `open_label_subjective` | Open-label design with a subjective / assessor-dependent primary endpoint | Assessment bias inflates observed effect |
| `non_ph_risk` (TTE) | Mechanism (e.g. IO) or precedent suggests delayed separation, but analysis is plain log-rank | Log-rank underpowered under non-proportional hazards; weighted tests or RMST may be needed |
| `ni_margin_wide` | NI margin larger than ~50% of the active control's historical effect | Margin may not preserve a clinically meaningful fraction of control benefit |
| `single_pivotal_thin` | Only pivotal study, modest N, p target ~0.05 | Unlikely to meet FDA's single-study substantial-evidence bar |
| `ssr_unblinded_no_firewall` | Unblinded sample size re-estimation without a stated firewall / DMC process | Operational bias risk; regulators will probe |
| `enrichment_unvalidated` | Powering relies on a biomarker-selected population but the biomarker assay / cutpoint isn't validated | Enrichment may not deliver the assumed effect amplification |
| `composite_driven_by_soft` | Composite primary where the powering effect is plausible only if the softest component drives it | Regulators and payers discount soft-component-driven composites |

## Disclosure / process flags

| Flag id | Fires when | Why it matters |
|---|---|---|
| `assumptions_undisclosed` | Sponsor states N and power but not the effect / control / SD assumptions | The implied delta is the only handle; note that the audit is by inference |
| `assumptions_changed_quietly` | Powering footnote differs across deck versions or vs the protocol paper | Signal of internal re-baselining |
| `enrollment_behind` | CT.gov enrollment or sponsor updates show pace materially behind plan | Calendar pressure leads to underpowered or delayed readouts |

When precedent is thin (rare disease, first-in-class), say so instead of
forcing a flag. "Unanchored powering" is the honest finding.
