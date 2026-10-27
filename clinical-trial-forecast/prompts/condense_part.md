READ FIRST — the style guide below overrides any conflicting instruction in this prompt.
{{style_guide}}

You are condensing one part of a {{report_kind}} on "{{indication}}" for {{asset}} to a HARD ceiling of {{max_words}} words including tables (it is now {{n_words}} words — you must remove at least {{must_cut}} words; an output that is not shorter, cuts too much, drops more than half the figures or does not start with a heading is rejected and the part keeps its current length; one still over the ceiling gets another round, up to three). Keep it all crisp and use bullets where possible. You have NO tools and must not add facts or change any number, parameter name, citation or figure reference you keep.

Where the forecast's uncertainty sits (allocate the words by load):
{{load_bearing}}

TERMINOLOGY: write for a reader with a pharma development / biostatistics modelling background. Standard terms of the trade need no explanation and should be used — responder difference in percentage points, non-responder imputation, estimand, co-primary, alpha allocation, attenuation, winner's curse, tetrachoric correlation, mixture prior, conditional power, Monte Carlo interval. What may not appear is the model's own private vocabulary: parameter ids, symbols (θ, K, W, R), labels the model gave a component or a scenario — say what the thing is and, for a scenario, what it assumes.

How to cut, in order:
1. Keep, untouched in substance: the headline numbers with their intervals; every table and figure that grounds a load-bearing step (the ones whose parameters top the ranking above); the Anchored on / Model choice / Why not / Swing bullets for those steps; the closing "What the model says" bullets (if present); every equation in an appendix.
2. Merge every stage or subsection that moves the estimate by less than ~2 pp into ONE short subsection of 3-5 bullets total ("Stages that are not drivers") — one sentence per stage with its value and its swing.
3. Delete: restatements of anything said earlier or in another part, background a specialist knows, throat-clearing, meta-commentary, duplicated tables (keep the most complete one), rows of a table that carry no weight (keep the load-bearing rows and say "n further trials, listed in the parameter page" if you drop rows), captions longer than one sentence, and any figure whose takeaway is already in another figure.
4. Compress prose into bullets: one full sentence per bullet (18-30 words), numbers with their conditions and citations kept — but never compress sentences into fragments or parenthetical number stacks, and never replace a plain-language name with a parameter id or symbol. Every list item goes on its own line with a blank line before the list (an inline " - item - item" run does not render).
5. A Summary part is condensed only by cutting whole bullets or sentences of lower load, never by shortening the sentences that remain.
6. If still over the ceiling, shorten the remaining subsections proportionally to their load, never the headline table.

Keep the original heading and heading levels; keep `![…](figures/…)` lines you retain exactly as they are; keep bracket citations. Return `section_markdown` only.

PART TITLE: {{title}}

INPUT:
{{section}}
