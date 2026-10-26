# Waypoint schemas

All waypoint files live in `<workdir>/waypoints/`. They are the resumable
state for the skill — read them on startup, write them after each phase.

## `progress.json`

Loop control (REQUIRED — write after every phase / revision).

```json
{
  "complete": false,
  "output_file": null,
  "current_phase": "enumeration",
  "iteration_notes": "What was accomplished this iteration"
}
```

When complete:

```json
{
  "complete": true,
  "output_file": "competitive_map.md",
  "current_phase": "synthesis",
  "iteration_notes": "What was accomplished this iteration"
}
```

## `enumeration.json` (Phase 1)

All competitors found.

```json
{
  "subject_asset": {
    "name": "<asset>",
    "company": "[company]",
    "indication": "[discovered or provided]",
    "target": "[primary target]",
    "mechanism": "[mechanism class]",
    "stage": "[development stage]",
    "key_data": "[headline clinical data]",
    "stated_timeline": "[company-stated milestones, or Not publicly available]",
    "sources": [{"source_url": "[URL]", "source_type": "[type]", "quote": "[verbatim]"}]
  },
  "competitors": [
    {
      "name": "[drug name]",
      "company": "[company]",
      "target": "[target]",
      "mechanism": "[mechanism]",
      "stage": "[Preclinical/Phase I/Phase II/Phase III/Under review/Marketed/Withdrawn/Discontinued]",
      "stage_source_url": "[URL]",
      "category": "[same_target|target_specific_soc|baseline_soc|adjacent_mechanism|other_pipeline|withdrawn_failed]",
      "threat_level": 1
    }
  ],
  "threat_summary": {
    "level_1": ["drug1", "drug2"],
    "level_2": ["drug3"],
    "level_3": ["drug4"],
    "level_4": ["drug5"]
  },
  "anomaly_flags": ["what failed or couldn't be found, and why it matters"]
}
```

## `findings.json` (Phase 2)

Accumulated research findings (written/updated in research phase).

```json
{
  "landscape": [
    {
      "topic": "[soc_outcomes|unmet_need|combinations|failure_reason]",
      "claim": "[finding]",
      "source_url": "[URL]",
      "source_type": "[type]",
      "quote": "[verbatim]"
    }
  ],
  "competitors": [
    {
      "name": "[drug name]",
      "aliases": ["[code name or brand name]"],
      "company": "[company]",
      "stage": "[stage]",
      "timeline": "[stated milestones and who's ahead of the subject asset]",
      "mechanism": "[detailed mechanism of action]",
      "efficacy": [
        {
          "primary_endpoint": "[endpoint]",
          "population": "[population and line of therapy]",
          "comparator": "[comparator]",
          "result": "[specific numbers with CI if available]",
          "trial": "[NCT# or trial name]",
          "head_to_head": "[true/false]",
          "source_url": "[URL]",
          "source_type": "[type]",
          "quote": "[verbatim]"
        }
      ],
      "safety": {
        "key_aes": [{"ae": "[name]", "incidence": "[X%]"}],
        "black_box": "[true/false/none]",
        "source_url": "[URL]",
        "source_type": "[type]",
        "quote": "[verbatim]"
      },
      "limitations": "[practical barriers, patient populations not served]",
      "positioning": "[line of therapy, typical patient]",
      "commercial": "[if marketed: publicly reported sales or market share]",
      "differentiation_vs_subject": "[how it compares to subject asset]",
      "other_sources": [
        {
          "supports": "[field the source backs: mechanism, limitations, commercial...]",
          "source_url": "[URL]",
          "source_type": "[type]",
          "quote": "[verbatim]"
        }
      ]
    }
  ],
  "coverage": {
    "level_1_researched": 5,
    "level_1_total": 5,
    "level_2_researched": 8,
    "level_2_total": 10,
    "gaps": ["drug X - missing efficacy data"]
  },
  "anomaly_flags": ["what failed or couldn't be found, and why it matters"]
}
```

For licensed records, follow the `licensed` exception in
`00-research-standards.md` in place of `source_url` and `quote`.

## `competitive_map.md` (Phase 3)

Final output. Markdown structure is specified in
`references/03-synthesis.md`.
