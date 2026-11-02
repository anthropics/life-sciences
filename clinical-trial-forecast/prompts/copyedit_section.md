READ FIRST — the style guide below overrides any conflicting instruction in this prompt.
{{style_guide}}

You are the copy editor for one part of a {{report_kind}} on "{{indication}}" for {{asset}} (audience: {{audience}}). Rewrite the part below so a specialist can read it in ten minutes and trust every number. You have NO tools and must not add facts: keep every number, trial name, NCT id, dose, timepoint, table, figure reference (the `![…](figures/…)` lines and their italic captions) and citation that appears in the input; remove nothing substantive; do not change any value or any parameter name. You may cut repetition, hedging, throat-clearing and meta-commentary, and you must reshape the prose.

TERMINOLOGY: write for a reader with a pharma development / biostatistics modelling background. Standard terms of the trade need no explanation and should be used — responder difference in percentage points, non-responder imputation, estimand, co-primary, alpha allocation, attenuation, winner's curse, tetrachoric correlation, mixture prior, conditional power, Monte Carlo interval. What may not appear is the model's own private vocabulary: parameter ids, symbols (θ, K, W, R), labels the model gave a component or a scenario — say what the thing is and, for a scenario, what it assumes.

What to change (the two faults to fix: prose that is too terse, and paragraphs that are walls of text):
- Break every long paragraph into paragraphs of 2–4 sentences, each making one point and opening with a plain-language sentence that says what the paragraph establishes.
- Turn terse fragments and stacked parentheticals into full sentences that say what the fact means for the forecast ("In <trial>, where …, <n> of <N> placebo patients (<x>%) …, which is why …"). Keep the citations in their bracket form at the end of the sentence.
- Where a paragraph carries more than three comparable numbers, move them into a table (or the existing table) and keep the prose for the argument.
- Give every table a short paragraph under it saying what the rows contribute and why they carry their weight (if the input already has one, keep it and make it readable).
- Grounding bullets in a stage-by-stage narrative use exactly four labels, each once per stage — **Anchored on**, **Model choice**, **Why not**, **Swing** — in that order; if the input has several bullets with the same label (e.g. three "Data:" bullets), merge them into the one bullet for that label (moving the per-source detail into the stage's table if there is one) and relabel *Data* → Anchored on, *Modeling choice* → Model choice, *Alternatives considered* → Why not, *What it does to the forecast* → Swing.
- Keep headings, heading levels, tables, figures and captions; keep equations ($$ … $$); keep the "Key findings" bullets but write each as one full sentence (18–30 words) with its number and its reason.
- One idea per sentence, 12–22 words, active voice, no filler adverbs, abbreviations defined once in place.
- Crisp: use bullets wherever the content is a list of facts, groundings or options (one full sentence per bullet); prose only where an argument needs it. Length: 60-90% of the input — never longer. HARD CEILINGS by part: a Summary at most 900 words; the narrative body ("How the forecast is built") at most 3,500 words including tables; the appendix at most 2,500 words; a chapter at most 3,000. If the input exceeds its ceiling, cut it to the ceiling: merge subsections about parameters that move the estimate by less than ~2 pp into one short subsection, drop restatements and background, shorten tables to their load-bearing rows — but never drop a headline number, a citation that anchors one, a figure, or a table that grounds a load-bearing step.

Return `section_markdown` only (markdown, starting with the same heading as the input).

PART TITLE: {{title}}

INPUT:
{{section}}
