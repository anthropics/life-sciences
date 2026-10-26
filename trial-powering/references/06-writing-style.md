# Writing Style

Write in the style of an experienced biostatistician or clinical development reviewer presenting the facts clearly and completely.

Open with the verdict stated plainly: the conclusion, the one number that
carries it, and the condition under which it flips. Then walk assumption by
assumption: name the assumption, show the worked number, state what changes
if it fails.

The expository model is Wittes 2002 (Epidemiol Rev 24:39-53). Write short
declarative sentences that state a number and say what it means, for example
(written for this skill, not quotations):

> In an event-driven trial, information comes from events, not from patients.

> Halving the assumed effect roughly quadruples the sample size.

> At 330 events the trial has 90% power for a hazard ratio of 0.70 and about
> 53% for 0.80.

## Citations

Cite every factual claim: quantitative data, clinical results, dates,
percentages, trial IDs, regulatory actions, safety data, quotes. Don't cite general knowledge or
your own interpretations.

Format: inline markdown links, [descriptive claim](url), placed on the
claim itself. Use numbered references only when a source appears 5+ times
or multiple sources support one claim.

Density: at least one citation per factual paragraph; more when claims in
the same paragraph come from different sources. Prose flow wins: group
adjacent claims from one source under a single link at the end of the
clause or sentence rather than linking mid-sentence fragments.

Link to specific data (trial pages, press releases, labels), not homepages,
and only URLs you actually fetched. Never fabricate statistics, dates,
names, or URLs; if a source fails or a number is not public, say so.

## Figures and tables

Tables for discrete values (audit rows, scenario grids, sample-size
options), units in headers. Charts only when shape carries the argument,
e.g. power vs true effect with the sponsor's assumption, MDE, break-even,
and precedent median marked on the curve. After rendering, Read the PNG and
re-render if any labels collide; an unreadable chart is worse than no chart.
