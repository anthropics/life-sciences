# Workflow Examples

These examples are adapted from the upstream `crispr-os` skill and reference set.

## Knockout Experiment

**User request:** `Design a CRISPR knockout experiment targeting TP53 in HEK293T cells`

Suggested workflow:

1. Read the biosafety guidance and confirm the target is appropriate.
2. Classify the request as a knockout workflow.
3. Review guide-design constraints for SpCas9 and NGG PAMs.
4. Run guide finding and ranking:

```bash
python3 scripts/pam_scanner.py <TP53_exon_seq> --cas SpCas9 --max-guides 5
python3 scripts/gc_content.py <guide_sequence>
python3 scripts/tm_calculator.py <guide_sequence>
python3 scripts/off_target_scorer.py <guide_sequence> <candidate_offtarget_sequence>
```

5. Generate a validation plan and provenance manifest:

```bash
python3 scripts/validation_planner.py --edit-type KO --cell-type HEK293T --gene TP53 --timeline-weeks 8 --json
python3 scripts/provenance.py --inputs <artifacts> --output manifest.json
```

## Base Editing

**User request:** `I need to convert an A-to-G at position 158 in BRCA1 exon 5`

Suggested workflow:

1. Classify the request as adenine base editing.
2. Find guides where the target adenine falls inside the editor window.
3. Penalize guides with likely bystander edits.
4. Prefer the guide with the cleanest edit window and lowest safety risk.

## CRISPRa

**User request:** `Activate VEGFA expression in iPSCs using CRISPRa`

Suggested workflow:

1. Classify the request as CRISPRa.
2. Review TSS-targeting constraints.
3. Search the promoter or TSS-proximal region for compatible PAM sites.
4. Check delivery constraints for iPSCs before finalizing the protocol.

## End-to-End Experiment Run

```bash
python3 scripts/experiment_runner.py run "Knock out TP53 in HEK293T cells"
python3 scripts/experiment_runner.py run "Knock out TP53 in HEK293T cells" --mode full
```

These commands come from the upstream project and represent the authoritative runtime entrypoint.
