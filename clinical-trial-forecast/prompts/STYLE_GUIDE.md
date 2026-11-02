# Report style guide

This guide is loaded into every report-writing prompt (the narrative, the appendix, the summary, condensing, copy-editing
and the claim-ledger edits) before the task instructions, and it overrides anything else in the prompt. Write like a
careful human analyst explaining the case to a colleague, not like a machine emitting compressed facts.

## Writing like a human
- Full sentences that say what a fact means. Never a terse fragment or a stacked parenthetical. Aim for 18–30 words
  per evidence sentence.
- Short paragraphs. A paragraph makes one point in 2–4 sentences and then stops. Never a wall of text: if a paragraph
  runs past five sentences or carries more than three numbers, split it, or move the numbers into a table.
- Each section and each paragraph opens with a plain-language sentence that says what it establishes and why it
  matters, before any number. Then the evidence, as complete, self-contained sentences with the connective tissue
  ("because", "so", "which means") that links one fact to the next. A reader should never meet a number without
  knowing what trial it comes from and why it is being brought up.
- One idea per sentence; prefer 12–22 words; active voice; name the actor ("the sponsor registered", "we estimate").
  Replace nominalisations with verbs. No filler adverbs ("notably", "critically", "importantly"), no rhetorical
  questions.
- Sentence limits: no sentence over 30 words, and a section's mean sentence under 22 words. No semicolon chains: a
  semicolon may join two clauses, never three. No parenthetical aside longer than eight words. A list of more than
  three items (countries, dates, trials) is a sub-bullet list or a table, never an inline run inside a sentence. One
  idea per sentence: the fact, then a second sentence for what it means. The rules against fragments never license
  fact-stacking: "full sentences" means short complete sentences, not long compound ones. These limits take precedence
  over the sentence-length targets above.
- State the issue the way it bites, not as an abstract category. Test for every risk: can the reader say in one
  sentence what would go wrong in the trial?
- Context for a specialist, not a primer. The readers are the audience the request names (by default, drug-development
  decision makers and clinical-development specialists): do not explain what an estimand or non-responder imputation is. What is needed is the orientation (name the trial and
  the measurement, say why it is relevant here), so every sentence stands on its own.
- No machine shorthand in prose: no `[derived: …]`, `[estimate]` or `[reported: …]` tags, no arrows or slashes
  doing the work of words, and no self-audit lines ("Verified: …", "Misread? … both used correctly"). Misread checks
  happen in the verification pass, and if something is wrong it is fixed or dropped, never annotated. Give provenance
  in words: "reported in the trial's primary paper", "derived from our census of recent placebo arms (Appendix A)",
  "our estimate". Citations stay in the short bracket forms ([NCT…], [PMID …], [DOI …], [URL …], [guidance: …],
  [registry: …]); the renderer turns them into numbered footnotes. Full details go in an appendix.
- Define each abbreviation once, in place, the first time it appears; no glossary. Use standard terms of art only, as
  a drug-development and clinical-development audience uses them (Phase 2 dose-ranging study, proof-of-concept study,
  first-in-human, pivotal trial, estimand, ITT); never invent labels such as "decision trial" or "go bar".

## Structure and skimmability
- The report is a short Summary, then a single narrative, then an appendix, followed only by the model's parameter
  table and a short Methods note. The narrative follows the model's stages and shows how each is grounded: the data
  and the modelling choice behind each adjustment, with tables where several sources feed a number and figures that
  earn their place. The appendix carries the exact mathematics, parameters, scenarios and validation. There are no separate analytical chapters, no separate analysis
  essay and no methods primer. Keep it all crisp and use bullets where possible. The whole report is about 6–9k words
  including tables.
- Each stage of the narrative grounds itself with exactly four bullets, one per label: Anchored on, Model choice, Why
  not, Swing. Never several bullets with the same label in a row: when several sources anchor a stage, they go into
  the stage's table. The effect waterfall (phase 2 delta → selection correction → haircuts → retention → true
  phase 3 effect → observed) opens the narrative, and the stages follow its bars.
- Most important first. Every section and subsection leads with the finding that most changes the decision, then the
  evidence, then the caveats. Cut background a specialist already knows, restatements of earlier tables and lists kept
  only for completeness: if a paragraph could go without losing a decision-relevant fact, it goes.
- Bottom line first, with enough detail to stand alone (about 20% more than a bare bullet): the recommendation or risk,
  the key numbers and conditions behind it, and the confidence. Key-findings bullets are full sentences, not
  number-stacks.
- Tables for 3+ comparable numbers; prose for the argument. Every table gets a short paragraph under it saying what
  each row contributes and why it carries the weight it does. Visible gridlines, and no decorative boxes. A table
  wider than about 7 columns is split into two stacked tables so it fits a page.
- The narrative is an argument told in prose, with tables and figures in support. It is not a table of numbers with
  connective fragments between them. Keep the evidence density (dataset tables, adjustments, arithmetic, citations),
  but carry it in short paragraphs of full sentences that explain what each number means for the forecast, with the
  tables doing the work of holding the numbers.
- Bullet lists: the label line is its own line; the list starts on the next line with a blank line before it; one
  full sentence per bullet. Lists render only when each item is on its own line, never as an inline " - a - b" run.
- Heading levels are consistent: report parts are level 2, each trial or stage subsection level 3; labels inside a
  block are bold run-in text, not headings.
- Every figure caption says in one sentence what to take from the figure; a figure without a takeaway is dropped. Use
  figures where they help (placebo response across trials, attrition, powering sensitivity, timeline).
- Short and skimmable: the reader sometimes reads on a phone. Density is not a virtue in itself. The goal is that a
  specialist can read a section in ten minutes and trust every number.

## Attention proportional to load
- The report's space and prominence follow what is load-bearing for the forecast's uncertainty, not the model's
  mechanics. The derivation of the forecast trial's effect from the earlier data is the centre of gravity of the
  Summary, the narrative and the appendix: the observed phase 2 delta and its sampling error, the selection
  (winner's-curse) correction, the phase 2 → phase 3 retention or attenuation and its calibration on precedent pairs,
  the design-change penalties, the translation to a second co-primary endpoint, and the reference arm for a
  head-to-head comparison, each where it applies. Lay out the evidence, the alternatives considered and the
  sensitivity of the estimate to each step in full.
- Powering at the expected effect (when conditional power is about 0.99), Monte Carlo error, the two-seed check and
  implementation housekeeping each get one sentence, or a line in a table. They are never a subsection of their own
  unless they actually move the estimate. A useful test: the words spent on a parameter should be roughly
  proportional to its swing in the tornado; anything below about 2 percentage points of swing is a table row, not a
  section.

## Language for an outside reader
- The Summary and the narrative prose may not use parameter ids, Greek symbols, mixture or component labels, scenario
  codes or any other name the model invented, because a reader outside the project cannot follow them. Name every
  quantity by what it is and explain every scenario in a clause. Ids and symbols belong in the appendix and the
  parameter table only, each introduced with its plain name once.
- Summary depth: the "What drives it" bullets are 2–4 full sentences each: what the quantity is, what evidence pins it
  (named trials with their numbers), and what it does to the answer (the range). No compressed clauses, no
  parenthetical number stacks, no " — " chains. Condensing must not buy length by turning sentences into fragments.
- Terminology: write for a reader with a pharma development or biostatistics background. Standard terms of the trade
  need no explanation and should be used: responder difference in percentage points, non-responder imputation,
  estimand, co-primary, alpha allocation, attenuation, winner's curse, mixture prior, conditional power, Monte Carlo
  interval, MMRM, retrieved-dropout imputation, tipping-point analysis. What may not appear is the model's own private
  vocabulary.

## Evidence and calibration
- Every headline claim is set in its evidence and in real-world examples: the primary source with its exact figure
  and condition, and 2–3 named trials or programs (with ids) that actually did this or showed this number. The
  reader's first questions are "is that real? have others done it? could this be a misread?". Answer the first two
  in the text; the third is handled by the verification pass.
- Conditions travel with numbers, never the number alone. Numbers carry unit and timepoint.
- Never launder a number through "project knowledge": carry the original primary citation forward.
- Calibration: no certainty words ("will", "cannot", "unreachable") on estimates or on medium- or low-confidence
  evidence. Confidence (H/M/L) and the load-bearing assumption go on every headline claim; give ranges, not points, for
  derived numbers; the text may not be more confident than the material it rests on.
- Facts versus ideas: the reader must always be able to tell a cited fact from our estimate, an assumption or the
  sponsor's plan. A registered design element or a press release's framing is the sponsor's disclosure, never an
  observed fact. Our own sequencing or timing must never read as if it were standard practice. How to signal it is a
  judgement call: "we propose", "we estimate", or framing a passage as a plan.

## The copy-edit pass
- The copy-edit pass may lengthen a part by up to 10% when that is what splitting sentences costs.
