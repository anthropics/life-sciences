---
name: igv-reports
description: Use when the user wants an HTML, clickable, browseable, offline, or emailable viewer of genomic data — phrases like "HTML IGV report", "offline IGV", "self-contained HTML", "clickable viewer", "create_report", "igv-reports", "email this viewer", or any browseable HTML of reads at variants, fusion breakpoints, SV junctions, viral integrations, ChIP peaks, ROIs, or ONT 5mC/5hmC methylation views at promoters/gene bodies/DMRs. Trigger even when the user doesn't say "igv-reports" — giveaway is HTML/clickable/offline plus genomic regions. Also fire on /igv-reports. DO NOT use for static PNG/PDF/SVG IGV screenshots — use the igv-screenshots skill instead.
---

# igv-reports

This skill builds **self-contained HTML genomic-region reports** with
[igv-reports](https://github.com/igvteam/igv-reports) (`create_report`).
Each report is a single browseable HTML containing the igv.js viewer plus
embedded data slices for every region. No server, no internet, no IGV
install needed at view time.

The skill has three entry points:
- **build** — one-shot: sites BED + BAM(s) ± VCF → HTML.
- **cohort** — multi-sample driver from a samplesheet → per-sample HTMLs + index.
- **prep-track** — utility: convert plain-gzip GFF/GTF/BED.gz into a
  bgzip + tabix-indexed track that igv-reports can load.

## What this skill is (and is not)

This skill is a **driver layer** on top of the upstream `igv-reports`
Python package by the IGV team
([github.com/igvteam/igv-reports](https://github.com/igvteam/igv-reports)).
The naming is unavoidable — both share the `igv-reports` name.

| Component | Source | Role |
|---|---|---|
| `create_report` CLI | upstream PyPI package `igv-reports` | does the actual HTML rendering |
| `scripts/build_igvreports.py` | **this skill** | wraps `create_report` with default-track resolution, cohort/samplesheet mode, SIF auto-detect |
| `scripts/verify_{report,cohort,anchors}.py` | **this skill** | post-render structural + content audits (not in upstream) |
| `scripts/prep_track.sh` | **this skill** | bgzip+tabix utility for annotation tracks |

## Install

```bash
# 1. Install the UPSTREAM igv-reports package (provides `create_report`):
pip install -U 'igv-reports>=1.16.0'

# 2. The skill's wrapper scripts ship inside this plugin. Once the plugin is
#    installed, the scripts live alongside SKILL.md.
```

If you only need raw `create_report` (no cohort mode, no verifiers, no
auto-tracks), skip this skill entirely and use upstream directly —
see [igvteam/igv-reports](https://github.com/igvteam/igv-reports) docs.

## Quickstart

```bash
python scripts/build_igvreports.py \
    --genome hg38 \
    --sites sites.hg38.bed \
    --bam tumor.bam normal.bam \
    --fasta /path/to/hg38.fa \
    --no-default-tracks \
    --extra-track /path/to/cpg_islands.bed.gz \
    --extra-track /path/to/gencode.v47.annotation.gff3.gz \
    --output report.hg38.html
```

If you run many reports across the same genome build, set up a databases YAML
once (schema in `references/databases_config_paths.md`) and point
`$IGV_REPORTS_DB_CONFIG` at it — then `--fasta` and `--no-default-tracks`
become optional.

## Environment overrides

All optional. Set per-shell or in a project `.env`:

| Var | Effect |
|---|---|
| `IGV_REPORTS_DB_CONFIG` | Path to a databases YAML resolving `--genome` to FASTA + default tracks (see `references/databases_config_paths.md`) |
| `IGV_REPORTS_SIF` | Path to an `igv-reports` apptainer SIF (offline / HPC use). Galaxy depot: `https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7cba7a3_0` |
| `SAMTOOLS_SIF_DEFAULT` | Path to a `samtools` SIF (verifier only) |
| `IGV_REPORTS_BIND` | Colon-separated bind paths for singularity. Unset = no binds. |

Driver flags `--fasta` and `--no-default-tracks` let you skip the databases
YAML entirely without setting any env var. `--no-apptainer` forces the PATH
`create_report` path even on a SLURM node. The hermetic `tests/unit/` suite
runs anywhere with `pytest` + Python ≥ 3.10.

## When to use which entry point

| User request | Entry point |
|---|---|
| "Make an HTML for these 5 SV breakpoints in tumor.bam" | **build** |
| "Give me one HTML per patient for the cohort integration calls" | **cohort** |
| "create_report fails with 'not BGZF' on this gencode" | **prep-track** |

## Defaults (locked in)

- Tracks always loaded, top-to-bottom in the viewer:
  1. CpG islands (BED, plain or bgzipped)
  2. Gencode full annotation (GFF3.gz, **transcripts + exons + CDS + UTRs**, NOT a gene-level-only file)
  3. RepeatMasker (BED.gz, bgzipped + tabix-indexed)
  Plus the user's BAM(s), VCF, and any extra tracks they pass.
- `--flanking 300` bp on either side of each site (good for SV breakpoints
  and point variants alike). Override per call if needed.
- `--standalone` so the HTML is offline-viewable.
- Output filename includes the genome tag — e.g. `cohort.hg38.html` —
  so downstream genome-tag enforcement hooks pass.
- Reference FASTA is resolved either explicitly via `--fasta` or via the
  YAML pointed to by `$IGV_REPORTS_DB_CONFIG`. Supported genome IDs out of
  the box: `hg38`, `mm10`, `mm39`, `t2t_CHM13v2_plusY`, `GRCh37` (extend
  `GENOME_ALIASES` in `scripts/build_igvreports.py` if needed).
- Per-genome default track availability when using the YAML is recorded in
  `references/databases_config_paths.md` — read it before assembling tracks
  so the skill doesn't try to load a track that doesn't exist for the
  selected genome.

## Sites BED format (critical)

igv-reports' BED parser reads fields **by position** and trips on a header
row (`ValueError: invalid literal for int() with base 10: 'start'`). Always
emit a **plain headerless 4-column BED**:

```
chr    start    end    name
chr2   25227855 25342590 DNMT3A_full_gene
```

Tab-separated. The `name` becomes the row label in the report's variant
table — make it specific enough to identify the site after deduping.

The project's `enforce-genome-tag.sh` hook requires a genome tag in the BED
filename: use `sites.hg38.bed`, not `sites.bed`.

## Pitfalls (the skill should encode and/or detect these)

| Symptom | Root cause | Fix |
|---|---|---|
| `ValueError: invalid literal for int()` on first row | Header row in sites BED | Strip header — plain BED |
| `UnicodeDecodeError: byte 0x8b` reading a track | igv-reports reading bgzip as text | Filename must end `.gff3.gz` / `.bed.gz` AND be true bgzip (check with `file <name>` for "extra field") |
| `tabix: not BGZF` | Track was plain-gzipped, not bgzipped | Run **prep-track** entry point |
| `tabix: out of order` while indexing | GFF/GTF/BED records not pos-sorted within chr | **prep-track** does `sort -k1,1 -k4,4n` before bgzip |
| Annotation track empty in viewer | Tabix returns no rows in displayed window — often correct biology (e.g., CGI-distal site). Confirm with `tabix file region` |
| Genome ID lookup fails with `--genome hg38` | igv.js bundled IDs require internet at view + render time. Use `--fasta /path/to/local.fa` instead (always works offline) |

Full pitfalls + create_report flag reference in `references/best_practices.md`.

## How to run — quick recipe

Ensure `create_report` is on PATH (`pip install -U 'igv-reports>=1.16.0'`).
If you use a conda env, activate it first.

Then call the bundled driver script (paths relative to the installed plugin):

```bash
python scripts/build_igvreports.py \
    --sites results/run/inputs/sites.hg38.bed \
    --bam tumor.bam normal.bam \
    --vcf calls.vcf \
    --genome hg38 \
    --fasta /path/to/hg38.fa \
    --no-default-tracks \
    --output results/run/reports/cohort.hg38.html
```

The driver:
- Resolves the genome's CpG / gencode / rmsk paths from `$IGV_REPORTS_DB_CONFIG`
  if set, skipping any that aren't configured for the chosen genome.
- Validates the sites BED is headerless and that all rows have `start < end`.
- Calls `create_report` with `--flanking 300 --standalone`.
- Writes a logs/ entry capturing the full command, the flanking value, the
  per-region embedded data sizes, and the resolved track list — useful for
  reproducibility and audit-trail expectations.

For multi-sample cohorts, use `--samplesheet samplesheet.tsv` instead of
`--bam/--vcf`. Samplesheet format: `sample, bam_tumor, bam_normal, vcf, sites_bed`.
The driver emits one HTML per sample plus a top-level `index.html` that lists
all samples with links. Layout matches the ATLL viral-integration reference
implementation:

```
results/<run>/
├── inputs/<sample>/sites.<genome>.bed
├── reports/<sample>.<genome>.html
├── reports/index.html
└── logs/run_<timestamp>.log
```

## prep-track — fixing a non-bgzip track

If a GFF3/GTF/BED.gz is plain-gzip rather than bgzip, igv-reports fails
silently or with an obscure error. Convert in place with backup:

```bash
bash scripts/prep_track.sh /path/to/track.gff3.gz
```

The script:
1. Backs up the original to `<name>.bak.original_gzip`.
2. `gunzip -c`s the file.
3. Sorts by `chr` then numeric `pos` (`sort -k1,1 -k4,4n`).
   (Gencode delivers records interleaved by feature type at the same locus —
   tabix requires pos-sorted.)
4. `bgzip`s in place.
5. `tabix -p <gff|gtf|bed>`s.
6. Verifies a sample tabix query returns rows.

Requires `bgzip` and `tabix` from htslib on PATH.

## When generating an answer.md / run.sh for the user

The driver script (`build_igvreports.py`) deliberately abstracts the
underlying `create_report` flags — it sets `--standalone`, `--fasta`, the
`--flanking 300` default, and the YAML-resolved annotation tracks
internally so the user doesn't have to remember them. That abstraction is
good for ergonomics but bad for auditability: a reviewer reading the
`answer.md` later can't see what flags are actually being invoked without
opening the driver source.

To keep both: when you produce a runnable command for the user, **also
include a code block titled "Equivalent direct create_report invocation"
that shows the fully-expanded command** with all flags and resolved track
paths inline. The user should see the wrapper command they're going to
run AND the underlying command it expands to. Example:

````
## Run

```bash
python build_igvreports.py --genome mm10 --sites peaks.mm10.bed \\
    --bam ./data/ip.bam ./data/input.bam \\
    --output reports/peaks_qc.mm10.html
```

### Equivalent direct create_report invocation

```bash
create_report peaks.mm10.bed \\
    --fasta /path/to/mm10.fa \\
    --flanking 300 --standalone \\
    --tracks ./data/ip.bam ./data/input.bam \\
        /path/to/mm10_CpGIslands.bed \\
        /path/to/gencode.vM25.annotation.gtf.gz \\
        /path/to/rmsk_all_repeats_mm10.bed.gz \\
    --title "ChIP-seq peak QC (mm10) — IP vs Input" \\
    --output reports/peaks_qc.mm10.html
```
````

This costs you ~10 lines and gives the reviewer a full audit trail. For
cohort runs, show the expanded form for ONE representative sample only —
the others differ only in BAM/VCF paths.

## Post-render verification

`scripts/verify_report.py` parses a built HTML and confirms it actually
contains what its inputs declared. Six checks: `html_exists`,
`html_min_size`, `region_count` (tableJson rows == sites BED rows),
`region_coords` (each BED row finds a matching `(chrom, start+1, end[, name])`
in tableJson — BED is 0-based, the HTML stores 1-based start), `region_sessions`
(sessionDictionary has one entry per row), and `tracks_present` (every
`name` from `--track-config` or every basename from positional `--tracks`
appears in the decoded igv.js session's `tracks[].name` list).

```bash
python scripts/verify_report.py \
    --html         results/<run>/reports/sample.hg38.html \
    --sites        results/<run>/inputs/sites.hg38.bed \
    --track-config results/<run>/inputs/tracks.json \
    --min-size-mb  1.0 \
    --out          results/<run>/reports/sample.verify.tsv \
    --fail-on-fail
```

Output is a TSV with columns `check / status / observed / expected / details`
(also printed to stdout). With `--fail-on-fail`, exits nonzero if any check
is FAIL — wire this into Snakemake / CI so the pipeline gates on render
quality, not just on `create_report`'s exit code.

NOTE: `--standalone` replaces every track URL with an inlined `data:` URL
after slicing, so URL paths are unrecoverable from the embedded session.
The check matches on track NAMES (which `--standalone` preserves) — for
`--track-config` JSON pass meaningful names; positional `--tracks` mode
uses basenames.

### Cohort-level verification (`verify_cohort.py`)

The per-sample verifier above confirms each HTML is internally consistent
but cannot tell whether sample-1's HTML accidentally embeds sample-2's BAM
(e.g., samplesheet typo, copy-paste, tumor/normal slot swap). For cohort
runs, `scripts/verify_cohort.py` adds five cross-sample checks:

| Check | What it asserts |
|---|---|
| `cohort_html_coverage` (global) | Each samplesheet row has exactly one HTML; flags missing + extras |
| `sample_tracks_match` (per-sample) | Each HTML's session contains every BAM/VCF basename declared in THAT row |
| `no_cross_sample_contamination` (per-sample) | Each HTML contains no basename that belongs to a DIFFERENT row's track columns (default tracks from `databases_config.yaml` are allow-listed) |
| `sample_id_embedded` (per-sample) | The `sample` column value appears in the HTML's `<title>` or filename |
| `index_consistency` (global) | `index.html` links exactly the samplesheet sample set; each target exists and is non-empty |

**Auto-invoked by default** at the end of `build_igvreports.py --samplesheet`
cohort runs. Disable with `--no-verify`; gate the pipeline with
`--fail-on-fail`. Standalone invocation:

```bash
python scripts/verify_cohort.py \
    --samplesheet samplesheet.tsv \
    --reports-dir results/<run>/reports/ \
    --genome hg38 \
    --out results/<run>/reports/cohort_verify.tsv \
    --summary results/<run>/reports/cohort_verify.summary.md \
    --fail-on-fail
```

The TSV adds a `sample` column on top of the per-sample verify schema, with
`"*"` for cohort-global rows. The markdown rollup (`--summary`) groups
PASS/FAIL counts by check + lists every failure inline.

Worked regression: `tests/integration/cohort_verify/scenarios.sh` builds a
3-sample cohort and asserts each of four corruption scenarios (missing
HTML, sample swap, index drift, truncated HTML) triggers the expected
check FAILs.

### Content verification (`verify_anchors.py`) — opt-in, slow

`verify_cohort.py` proves the HTML *says* the right thing. It can NOT
confirm the embedded BAM *slice* contains the data it claims to. Two
failure modes slip past structural checks:

1. **Sample swap with matching basename** — the cohort loop wired the wrong
   BAM into `sample_1`'s build, but the swapped BAM's `Path.stem` happens
   to match what `sample_1`'s row declared (or two files in different dirs
   share a basename). Track name passes; slice content is wrong.
2. **Silent empty slice** — region rendered, but the slice has 0 reads
   (failed `samtools index`, source BAM corruption, coords outside coverage).

`scripts/verify_anchors.py` closes the gap by re-running `samtools view -c`
against both the source BAM (at generate time) and the embedded slice (at
verify time), then comparing counts. Two-mode workflow:

```bash
# 1. After the cohort renders cleanly, freeze the read counts as a regression fixture.
python scripts/verify_anchors.py generate \
    --samplesheet samplesheet.tsv \
    --sites sites.hg38.bed \
    --out anchors.hg38.tsv

# 2. Re-verify any time after — works against a fresh build of the same inputs,
#    or to audit an existing HTML for unexpected content drift.
python scripts/verify_anchors.py verify-cohort \
    --samplesheet samplesheet.tsv \
    --reports-dir results/<run>/reports/ \
    --genome hg38 \
    --anchors anchors.hg38.tsv \
    --out results/<run>/reports/cohort_verify_anchors.tsv \
    --fail-on-fail
```

Or chained into the build driver:

```bash
# Freeze anchors at build time:
python scripts/build_igvreports.py --samplesheet ... --anchors-mode generate \
    --anchors anchors.hg38.tsv

# Verify a later build against frozen anchors:
python scripts/build_igvreports.py --samplesheet ... --anchors-mode verify \
    --anchors anchors.hg38.tsv --fail-on-fail
```

Anchors TSV schema (`#`-prefixed header per lab BED convention):

```
#sample	track_name	chrom	start	end	expected	tolerance	min	max	notes
```

`tolerance` is a ratio (default 5%). `min`/`max` are absolute bounds that
override tolerance when set — useful for known-positive sites like
"this integration must have ≥20 reads".

samtools is resolved in this order: `--samtools-sif PATH` → `$SAMTOOLS_SIF`
→ `$SAMTOOLS_SIF_DEFAULT` → PATH `samtools`. On HPC, prefer a SIF to avoid
the NFS conda cold-start tax.

**Why opt-in and not default:** the verify step shells out to samtools per
(sample × region) and indexes each slice — ~1 s/anchor. For a 6-sample
cohort × 50 regions that's ~5 min on top of the structural verify (which
runs in seconds). Reach for this when sample swap or content regression
is a real concern; the structural verifier is sufficient for routine builds.

Worked regression: `tests/integration/anchor_verify/scenarios.sh` builds a
2-sample cohort and asserts each of four content scenarios (tolerance
violation, min-bound violation, corrupted slice, missing anchor) triggers
the expected PASS / FAIL / SKIP outcome.

## Output and workflow logging

Every run logs to `logs/run_<YYYYMMDD_HHMMSS>.log` next to the reports dir.
The log captures:
- Resolved track paths (per genome, after databases_config.yaml lookup).
- The exact `create_report` command.
- The flanking value used (default **300 bp** — this is the value that's
  baked into all the embedded data slices, so audit trails depend on it).
- Per-region embedded data sizes (extracted post-render so the user can
  see which regions inflated the HTML).
- Total HTML size.

This satisfies CLAUDE.md §"Logging and Audit Trail" — every run is
reproducible from the log alone.

## Track choice nuances

For gencode on hg38, the default points at
`gencode.v47.annotation.gff3.gz` (full annotation, bgzip + tabix). This
gives transcript models with exons / CDS / UTRs. The gene-level-only
companion (`gencode.v47.genes.annotation.sorted.gff3.gz`) renders only
solid gene boxes and is fine for high-zoom views, but the full annotation
is the right default for read-level inspection at integration / fusion /
SV junctions.

For mouse genomes, `databases_config.yaml` ships `.gtf.gz` paths instead.
GTFs work in igv-reports if bgzip + tabix-indexed; **prep-track** converts
plain-gzip GTFs the same way it does GFF3s.

For T2T-CHM13, only the FASTA + GTF + CGI are indexed in our DB; rmsk is
absent and is auto-skipped by the driver. The variant table will load
without rmsk; flag this in the run log.

## Common-case examples

The `examples/` directory has runnable templates:

- `single_sample.sh` — one BAM + one VCF + a sites BED → one HTML.
- `cohort_samplesheet.sh` — TSV-driven multi-sample run.
- `prep_track_demo.sh` — convert a plain-gzip gencode to bgzip+tabix.
- `methylation_ont/` — ONT 5mC/5hmC viewer (BAM with `colorBy: basemod2`
  + per-sample bedGraph at fixed y-axis 0..100). End-to-end worked
  example with pre-sliced data; recipe.md explains the slots.

These are reference implementations; copy and edit them for new runs
rather than starting from scratch.

## Tests

Three-layer suite under `tests/`, orchestrated by `tests/run_all.sh`:

| Layer | What it covers | Runtime | Needs |
|---|---|---|---|
| **unit** (`tests/unit/`) | parser layer of `verify_report.py` + `verify_anchors.py` — TSV loading, status decision, session-entry locator, balanced-brace JSON extractor, decode round-trip — all with synthetic inputs | ~1 s | pytest |
| **smoke** (`tests/smoke/`) | `samtools_count` / `samtools_index` / full slice-decode-and-count round-trip against the committed `tests/fixtures/tiny_colo829.hg38.bam` (457 KB, sliced from public ONT COLO829 release) | ~3 s | pytest + samtools (SIF or PATH) |
| **integration** (`tests/integration/`) | end-to-end: build a 2-/3-sample cohort, structural verify, anchor verify, run 4 corruption scenarios per verifier | ~7 min cold, ~30 s cached | full cohort BAMs (lab default OR `IGV_REPORTS_TEST_BAM_{1,2,3}` env override). SKIPs with exit 77 if neither is available |

```bash
bash tests/run_all.sh                  # all three layers
bash tests/run_all.sh --unit-only      # ~1 s — fastest feedback loop
bash tests/run_all.sh --no-integration # ~12 s — works on any machine
bash tests/run_all.sh --integration-only
```

The fixture provenance + regeneration recipe live in
[tests/fixtures/README.md](tests/fixtures/README.md). Anchor counts the
smoke layer expects (chr2=5, chr7=9) are the contract — any fixture
regeneration that changes them must also update the smoke test constants.

## ONT methylation viewers (specialized path)

For per-read 5mC/5hmC visualization the positional `--tracks` API does
not work — you need named tracks with `colorBy: "basemod2"` on the BAMs
and `min: 0, max: 100` on the bedGraph tracks (cross-sample y-axis lock,
see `rules/igv.md`). Use the `--track-config <json>` passthrough:

```bash
# 1. Write a YAML spec listing samples (see tracks_spec.example.yaml).
# 2. Generate tracks.json with the right defaults baked in:
python scripts/generate_tracks_json.py \
    --spec tracks_spec.yaml --run-dir results/<run>/ \
    --out results/<run>/tracks.json

# 3. Build the report:
python scripts/build_igvreports.py \
    --sites results/<run>/sites.hg38.bed \
    --track-config results/<run>/tracks.json \
    --genome hg38 --flanking 0 \
    --type mutation --info-columns name \
    --output results/<run>/methylation_report.hg38.html
```

Key methylation-specific defaults:
- `--flanking 0` (sites BED already encodes the window — promoter/gene span).
- `--info-columns name` (surface the BED `name` column in the variant table).
- `--type mutation` (one-locus view per row; not split-screen).
- bedGraph not bigwig — igv-reports cannot slice `.bw` directly.

When `--track-config` is set the driver bypasses the auto-resolved
default annotation tracks (CGI / gencode / rmsk) and the `--bam` /
`--vcf` / `--extra-track` flags — the JSON is the source of truth.
Build annotation slices into the JSON instead.

**`--apptainer` is auto-detected**: the driver flips to the apptainer SIF
pointed to by `$IGV_REPORTS_SIF` (igv-reports 1.16.0, ~83 MB, pulled from
the Galaxy depot) when `SLURM_JOB_ID` is in the environment — i.e. running
on a compute node where the NFS conda cold-start tax matters. On the login
node or when `$IGV_REPORTS_SIF` is unset, the driver uses PATH
`create_report`. Override with `--apptainer` / `--no-apptainer`; the
decision lands in the run log.

Full recipe and rationale: `references/methylation_ont.md`. Worked
example with real data: `examples/methylation_ont/`.

## See also

- `references/best_practices.md` — full create_report flag reference,
  format gotchas, performance notes. Read this if a run fails in a way
  not listed in the Pitfalls table above.
- `references/databases_config_paths.md` — per-genome track availability
  matrix and exact YAML keys. Read this when adding a new genome or
  diagnosing a missing-track warning.
- `references/methylation_ont.md` — ONT 5mC/5hmC cheat-sheet (colorBy,
  min:0/max:100, flanking=0, bedGraph vs bigwig, EPDnew lookup).
- `scripts/build_igvreports.py` — the driver. Reads `--samplesheet` or
  `--bam/--vcf` direct-args, resolves tracks, validates the sites BED,
  writes the HTMLs and the run log. Supports `--track-config <json>`
  passthrough for fully-styled track sets.
- `scripts/generate_tracks_json.py` — YAML spec → tracks.json with
  ONT-methylation defaults baked in (colorBy=basemod2, min:0/max:100,
  group-paired Okabe-Ito colors).
- `scripts/verify_report.py` — post-render structural verifier; parses
  the HTML's embedded tableJson + sessionDictionary, confirms region
  count / coordinates / track names match the inputs. Emits a verify.tsv
  and gates on `--fail-on-fail`.
- `scripts/verify_cohort.py` — cohort-level verifier; layered on top of
  verify_report's per-sample checks, adds cross-sample contamination
  scanning + index.html / sample-id consistency. Auto-invoked at the end
  of `build_igvreports.py --samplesheet`; standalone-runnable too.
- `scripts/verify_anchors.py` — content verifier; samtools-counts the
  embedded BAM slices and compares to anchors frozen from the source BAMs
  at build time. Catches sample swaps that share basenames and silent
  empty slices. Opt-in via `--anchors-mode generate|verify` on the build
  driver; slow (~1 s/anchor). See SKILL.md content-verification section.
- `scripts/prep_track.sh` — gunzip → sort → bgzip → tabix utility.
- `igv-screenshots` skill — the **static PNG/PDF/SVG** counterpart based
  on igver. Use it instead of this one when the deliverable is a
  publication-quality figure rather than a clickable viewer.
- Upstream development: https://github.com/sahuno/igv-reports-skill
  — file issues there for skill-level bugs; file issues at
  https://github.com/igvteam/igv-reports for `create_report` rendering bugs.
