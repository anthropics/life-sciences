You are setting up a probability-of-technical-success forecast for one registered clinical trial. This step replaces the project configuration a person would otherwise write by hand: everything later steps need to know about the trial, the indication, the report and the model comes from what you write here.

Read:
- question.md: the person's request, verbatim. It names the trial to forecast and may say what counts as success, which decision the forecast informs, or where to look.
- raw/ctgov/{{focal_nct}}.json: that trial's registry record.
- intake/notes/background.md, if it exists: background the person supplied. Treat it as confidential unless it is plainly public (a press release, a publication).

This step has no web access on purpose. Use only these files and what a specialist already knows; later steps widen the evidence from the run's files (and from outside sources, when the person allowed them). Do not narrow early.

Write project.json with these fields:
- `title`: a short report title naming the trial (acronym or NCT id), the drug and the indication.
- `scope_label`: one line naming the drug's mechanism class and the indication.
- `public_thesis`: the request restated for the web-enabled researchers, who never see question.md or the background: the trial, what is being forecast, and any guidance on where to look. Leave out anything the background marks as confidential.
- `focal_trial`: `nct_id` ({{focal_nct}}) and `name` (its acronym or short title).
- `indication`: `name`, `aliases`, `reference_population` (the population the trial enrolls, in a line), and `key_endpoint`: `label`, `description` (the efficacy measure this trial and its precedents are judged on), `keywords` (5-10 words or phrases that find that endpoint in registry outcome text) and `typical_timepoint_weeks` (a number, or null).
- `target`: the drug's biological target, `name` and `aliases`.
- `asset`: the drug under test. `public` is true for a drug named in a public registry record, which is the usual case; set it false only if the person asks for the drug to be kept confidential. `name`, `code` (a sponsor code, or ""), `developer`, `modality`, `stage` (for example "phase 3").
- `report`: `goal` (what the forecast is for, in a sentence), `decision` (the decision it informs: the person's, if stated; otherwise what an investor or sponsor would decide on it), `audience` ("drug-development decision makers and clinical-development specialists" unless the request names another), `scope_note` (what the report covers and leaves out, in a line), `flip_criteria` (what evidence would change the forecast materially), and `sections`: 4-7 analytical sections, each with `id` (snake_case), `title` and `question`. The sections are the questions the model's parameters rest on, so the research loop can fill them: typically the control-arm response under this design and population, how much of the earlier-phase effect will carry into this trial (the precedent pairs in the class and indication), the design, analysis and powering of this trial, safety, tolerability and operational risks that could end it, and the bar for success (regulatory and clinical). Adapt them to this trial; do not add a section the model would not use.
- `simulation`: `target_description` (the trial to simulate, in a line: design, arms, N, primary endpoint and timepoint, from the record), `success_definition` (what counts as technical success: the person's definition, verbatim, if question.md gives one; otherwise "the trial meets its primary endpoint(s) per the pre-specified analysis"), and `considerations` (3-8 specific things the parameter sheet must weigh for this trial, for example a change of population or dose between phases, or a co-primary endpoint).
- `registry`: the searches that build the evidence base.
  - `indications`: the trial's indication first, then up to 2 others only if they hold the precedent trials this forecast needs (for example a neighbouring indication where the same class ran phase 2 and phase 3). Each: `id` (snake_case), `name`, `aliases`, `ctgov_condition_query` (a ClinicalTrials.gov condition query, specific enough that a pull stays in the hundreds), `key_endpoint_description`, `key_endpoint_keywords`, `reference_population`, `why`.
  - `mechanism_classes`: 10-30 classes covering the drugs tested in these indications, at the granularity the literature compares drugs (the drug's own class and its competitors; an antibody and a small molecule against one target are separate classes). Do not add "other", "non-drug" or "unknown"; the tools add them.
  - `analog_features`: 3-6 features that make another program a precedent for this one regardless of indication.
  - `pathway_interventions`: 5-25 registry intervention names (INNs or sponsor codes, one drug per entry, spelled as the registry lists them; never a brand name, class name or dose) of drugs on the same target or pathway and the direct competitors. Their trials are pulled in any condition, because phase 2 to phase 3 pairs in the class are the model's anchor. Include the drug under test itself when it is public. Leave out old generic drugs with huge unrelated trial footprints.
- `notes`: how you read the request, in a short paragraph.

Every field except `notes` reaches the steps that search the web, so keep out of them anything the background marks as confidential.

Return JSON matching the project schema exactly.
