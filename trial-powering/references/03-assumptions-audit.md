# Phase 3: Assumptions audit

The power calculation is only as good as its inputs. This phase finds external
evidence for what each input *should* be, so Phase 4 can stress the design
against realistic values instead of arbitrary perturbations.

Every audited value must carry a citation (`source_url`, `source_type`,
`quote`) per `00-research-standards.md`. The output of this phase is the
evidence base that makes the sensitivity grid defensible rather than
hypothetical.

## The four assumptions that decide most trials

### 1. Treatment effect size

What the sponsor assumed vs what the evidence supports.

Sources, in rough priority:
- This asset's own Phase 1b/2 data in this indication (sponsor materials, conference abstracts, publications). Apply a **Phase 2 to Phase 3 shrinkage** expectation: open-label, single-arm, or biomarker-enriched Phase 2 effects routinely shrink 20-50% in randomized Phase 3. Say so explicitly.
- Same-target / same-MOA assets in the same indication (the class effect). The `precedent-analysis` skill returns these. Otherwise use the `pubmed` MCP and `WebSearch` for "[target] [indication] phase 3 [endpoint]".
- Meta-analyses of the class if they exist (`pubmed` MCP, `biorxiv` MCP for medRxiv preprints).

Record as: `{sponsor_assumed, asset_prior_data, class_median, class_range, sources[]}`.

### 2. Control arm performance

This is where the most Phase 3 misses happen. Standard of care improves,
trial populations are healthier than historical cohorts, and placebo response
inflates in CNS / pain / derm. A control arm that outperforms the assumption
shrinks the observable delta even if the drug works exactly as hoped.

Sources:
- Contemporary Phase 3 trials in the same line/population (last ~5 years). Pull the **observed control arm** result, not the assumption. `clinical_trials.analyze_endpoints` + PubMed.
- Disease registries / natural history studies if the comparator is placebo / BSC.
- For oncology TTE: recent control-arm medians in the same line, histology, biomarker subset.
- For CNS/pain/derm continuous endpoints: recent placebo-arm change-from-baseline values (placebo response has trended up over time in many indications).

Record as: `{sponsor_assumed, precedent_values: [{trial, year, value, source}], precedent_median, precedent_range}`.

### 3. Variability / overdispersion

For continuous: SD of the primary endpoint. For counts: overdispersion k.

Sources: published trials using the same instrument (e.g. ADAS-Cog, EDSS,
HAM-D, EASI) in a comparable population. SD is usually in Table 1 or the
primary results table. Sponsors often take SD from their own small Phase 2;
larger trials almost always show higher SD.

### 4. Dropout / event accrual

Sources: precedent trials' CONSORT diagrams (n randomized vs n analyzed),
sponsor's own enrollment updates (CT.gov history, press releases), and for
TTE the event-count updates if disclosed.

## Subagent fan-out

If subagents are available, run these three in parallel and collate (use a
domain research agent type if one is configured; otherwise general-purpose):

**Brief A (effect-size precedent):**
> For [asset / MOA] in [indication, line, population], find the observed
> [endpoint metric] from this asset's prior trials and from same-class assets.
> Return a table: trial, phase, N, observed effect, CI, source URL. Note
> whether each was randomized or single-arm. For each precedent trial, also
> locate the sponsor's deck covering its design or readout (Events &
> Presentations page; data-event decks preferred) and return the
> deck URL, date, and what it discloses; record not_found if it does not
> exist.

**Brief B (control-arm precedent):**
> For Phase 2/3 trials in [indication, line, population] reading out in the
> last 5-7 years, extract the control/placebo arm result on [endpoint]. Return
> a table: trial, year, control arm regimen, control arm result, N, source.

**Brief C (variability + dropout precedent):**
> For trials measuring [instrument / endpoint] in [population], extract the
> reported SD (or overdispersion) and the dropout rate. Return a table.

Collate into `<workdir>/assumptions_audit.json`:

```json
{
  "effect_size": {"sponsor": ..., "asset_prior": ..., "class_median": ..., "class_range": [..], "verdict": "in-line|optimistic|conservative", "sources": [...]},
  "control_arm": {"sponsor": ..., "precedent_median": ..., "precedent_range": [..], "verdict": "...", "sources": [...]},
  "variability": {"sponsor": ..., "precedent_median": ..., "verdict": "...", "sources": [...]},
  "dropout":    {"sponsor": ..., "precedent_median": ..., "verdict": "...", "sources": [...]},
  "other": [{"param": "...", "note": "..."}]
}
```

`verdict` is a plain-language call: `in-line` (within precedent IQR),
`optimistic` (sponsor assumption is more favorable than precedent median by a
material amount), or `conservative`. Reserve `optimistic`/`conservative` for
differences large enough to move power by >5 absolute points; otherwise call
it `in-line` and note the gap.

## When precedent is thin

Rare disease, first-in-class, novel endpoint: there may be no good precedent.
Say so. Use the sponsor's own earlier-phase data with explicit shrinkage, and
widen the sensitivity grid in Phase 4 to compensate. Thin precedent is itself
a finding (it means the powering is unanchored, which is a risk regardless of
perspective).
