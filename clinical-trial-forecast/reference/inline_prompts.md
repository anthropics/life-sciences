# The claim-extraction prompt

The claim-extraction step follows the prompt below. It reads the report from a file, and each claim also gets a public
restatement (the second paragraph below), because this skill's claim checkers may use outside sources.

## Claim extraction (the claims step; MAX_CLAIMS is the number your brief gives)

Extract the headline claims from this probability-of-technical-success forecast that a skeptical reader would want
checked: every bottom-line / summary / thesis bullet, every issue heading that asserts a fact, every quantitative fact
that carries a recommendation or an estimate (trial design facts, prior-phase results, placebo rates and attrition,
precedent effect sizes, attenuation factors, failure-mode rates, regulatory rules and timelines, costs), every stated
simulation parameter value with its justification (if the report has a model), and every recommendation. One entry per
checkable claim, verbatim text in `claim`, the surrounding 1-2 sentences in `context`, the section name, and the kind.
Up to MAX_CLAIMS claims, most decision-relevant first. Do not paraphrase.

The checkers work from public sources and never see the report, so for each claim also write `public_claim` and
`public_context`: the same checkable assertions restated so that they stand alone and nothing confidential leaves the
project (nothing that comes only from the person's request in question.md or the background they supplied in
intake/: a bar, a threshold, a definition of success or any other figure the request set is confidential). The trial,
the drug and every public fact keep their names and numbers exactly. If the claim's substance rests only on that
request or background, or is purely our own judgement with nothing public to check, set `sponsor_only` true.

REPORT: the whole of report/claims_source.md (the report up to its first Appendix heading).
