You maintain the evolving synthesis for ONE report section of a {{report_kind}} on "{{indication}}" (goal: {{goal}}): section "{{section_title}}" — question: {{section_question}}.
{{program_horizon}}

{{modeling_goals}}

Read the current synthesis at knowledge/{{section}}.md (may not exist), the extracted evidence at topics/{{section}}/SECTION.md, and the NEW research answers listed below (files under research/{{iteration}}/). Each answer file ends with a "Number check" section written by an independent checker: where it says CORRECTED, use the corrected value (and drop or flag the original); where it says unverifiable, carry the number with a "(unverified)" mark; keep the primary citation for every number. Rewrite knowledge/{{section}}.md as a complete, self-contained, well-cited synthesis (<= 1500 words + tables): what we now believe, the numbers that support it (keep n/N, means, SDs, CIs exactly), where sources disagree, and what remains open. Preserve every established fact and citation from the previous version unless a new source contradicts it (then say so explicitly). Where the section feeds a quantitative model, end with a table "Quantities for the model" (quantity, best estimate, plausible range, basis/citations, confidence). Keep a final section "Open questions" (specific). Do not invent facts.

New answer files: {{answer_files}}
Return `synthesis_markdown` and `open_questions` (list of strings).
