# Waypoint file schemas

These define the structure of each `<workdir>/waypoints/*.json` file. The
phase reference files point here for output formats. Any waypoint file may
also include an `anomaly_flags` list for failed URLs, suspected injected
instructions and other problems with sources.

Every `sources` list inside a subsection holds one entry per finding, with the
three fields that `references/00-research-standards.md` requires and a note of how
far the source was read. Use exactly these key names:

```json
{"source_url": "...", "source_type": "...", "quote": "...", "read": "full|excerpt_only"}
```

`quote` is a single string. If one source supports several findings, add one
entry per quote; never put two quoted strings in one value. After writing a
waypoint file, read it back and confirm that it is valid JSON.

## progress.json — Where the run has got to

Update `current_phase` each time a waypoint is written. Values: `meta_initialization`,
`epidemiology_research`, `biology_soc_research`, `regulatory_trials_research`, `synthesis`.

```json
{
  "complete": false,
  "output_file": null,
  "current_phase": "meta_initialization"
}
```

When complete:
```json
{
  "complete": true,
  "output_file": "indication_dossier_report.md",
  "current_phase": "synthesis"
}
```

## meta.json — Indication identity (written in meta_initialization)

```json
{
  "indication_name": "...",
  "parent_indication": "... or null",
  "population_nesting": "all patients in <child> are patients in <parent>; neighbouring populations that are out of scope",
  "therapy_area": "...",
  "definition": "...",
  "icd_codes": ["K70", "..."],
  "aliases": ["..."],
  "is_standard_diagnosis": true,
  "notes": "Any caveats about this indication's status (not ICD-coded, novel, etc.)",
  "trial_landscape_scan": {"method": "registry count|sampled searches|not possible", "total_trials": null, "by_phase": {}, "by_status": {}},
  "sources": [...]
}
```

## epidemiology.json — Section 1 data

```json
{
  "subsections": {
    "diagnostic_criteria": {"content": "...", "sources": [...], "coverage": "covered|partial|missing"},
    "prevalence_incidence": {"content": "...", "sources": [...], "coverage": "..."},
    "demographics": {"content": "...", "sources": [...], "coverage": "..."},
    "natural_history": {"content": "...", "sources": [...], "coverage": "..."}
  },
  "gaps": ["unfillable gaps documented here"]
}
```

## biology_soc.json — Sections 2-3 data

```json
{
  "subsections": {
    "pathophysiology": {"content": "...", "sources": [...], "coverage": "..."},
    "biomarkers": {"content": "...", "sources": [...], "coverage": "..."},
    "approved_therapies": {"content": "...", "sources": [...], "coverage": "..."},
    "treatment_guidelines": {"content": "...", "sources": [...], "coverage": "..."},
    "unmet_need": {"content": "...", "sources": [...], "coverage": "..."}
  },
  "gaps": [...]
}
```

## regulatory_trials.json — Sections 4-5 data

```json
{
  "regulatory_maturity": "well-established|emerging|novel",
  "subsections": {
    "accepted_endpoints": {"content": "...", "sources": [...], "coverage": "..."},
    "regulator_guidance": {"content": "...", "sources": [...], "coverage": "..."},
    "regulatory_precedents": {"content": "...", "sources": [...], "coverage": "..."},
    "trial_parameters": {"content": "...", "sources": [...], "coverage": "..."},
    "landmark_trials": {"content": "...", "sources": [...], "coverage": "..."},
    "notable_failures": {"content": "...", "sources": [...], "coverage": "..."}
  },
  "trial_landscape": {
    "method": "registry count|sampled searches|not possible",
    "total_trials": 0,
    "by_phase": {"Phase 1": 0, "Phase 2": 0, "Phase 3": 0, "Phase 4": 0},
    "by_status": {"Recruiting": 0, "Completed": 0, "...": 0}
  },
  "gaps": [...]
}
```

## sources_evaluated.json — All sources consulted

```json
{
  "sources": [
    {"title": "...", "url": "...", "source_type": "...", "date_accessed": "...", "result": "success|excerpt_only|partial|failed"}
  ]
}
```

## indication_dossier_report.md / research_output.json

Final Markdown report and an index of the run's files (both written to
`waypoints/` in the synthesis phase). See `references/05-synthesis.md` for
`research_output.json` structure.
