---
name: methylation_ont
genome: hg38 | mm10 | mm39 | t2t
assay: ONT 5mC + 5hmC (CpG)
worked_example: ../examples/methylation_ont/
---

# ONT methylation viewer — cheat-sheet

Targeted reference for building an igv-reports HTML that shows per-read
5mC/5hmC base-modification calls (BAM, basemod2 coloring) plus per-sample
methylation-fraction bedGraph tracks at fixed promoter / gene / DMR windows.

When this skill needs to build a methylation viewer, the **default path
(positional `--tracks`) is wrong** — methylation viewers need named,
colored, y-axis-locked tracks. The right path is:

```bash
build_igvreports.py --track-config tracks.json ...
```

with `tracks.json` either generated from a YAML spec (see worked example)
or hand-written from `tracks.template.json`.

## The four-thing checklist

### 1. BAM tracks need `colorBy: "basemod2"`

```json
{
  "name": "<sample>",
  "url": "<bam>",
  "indexURL": "<bam>.bai",
  "format": "bam",
  "type": "alignment",
  "colorBy": "basemod2",
  "showSoftClips": false,
  "displayMode": "COLLAPSED"
}
```

Without `colorBy: "basemod2"`, the BAM renders as plain alignments
without the per-base 5mC/5hmC colors that are the whole point of the
view. `displayMode: "COLLAPSED"` keeps the BAM panel short so the
bedGraph summary tracks below stay visible.

### 2. bedGraph tracks need fixed `min: 0, max: 100`

```json
{
  "name": "<sample> 5mC",
  "url": "<bedgraph>",
  "format": "bedgraph",
  "type": "wig",
  "color": "rgb(0,68,136)",
  "min": 0, "max": 100
}
```

modkit's bedmethyl output is **percent (0..100)**, not fraction (0..1) —
the y-axis ceiling must be 100. IGV's per-track autoscale defaults
differ per track and hide real cross-sample differences (one sample
might autoscale to 0..82, the next to 0..100; same bar height means
different methylation). Lock all samples' bedGraph tracks to the same
0..100 range. See `rules/igv.md` for the original incident.

**Use bedGraph, not bigwig.** igv-reports' Python slicer (`utils.getreader`)
dispatches on file extension and has no `.bw` reader — runs fail with
`Exception: Unknown file format`. Pre-slice bigwigs over the report
regions with `bigWigToBedGraph -chrom -start -end <bw> <bg>`, one
output per region, then `cat >>` them into a single bedGraph (UCSC
`bigWigToBedGraph` opens `/dev/stdout` with `O_TRUNC` between calls —
piping multiple invocations loses everything but the last region).

### 3. `--flanking 0` when sites encode the desired window

For methylation viewers the sites BED almost always carries the desired
window directly (a promoter span, a DMR, a gene body). Adding 300 bp of
flanking adds nothing and shifts the initial viewer frame. Pass
`--flanking 0` and let the BED row coordinates be the frame.

The 300 bp default is right for the SV/integration breakpoint workflow
this skill was extracted from — there the BED row is a one-base
breakpoint and you need flanking to see read support.

### 4. Sites BED with `#chrom\tstart\tend\tname` comment header is fine

The skill's older docs say "headerless" because non-`#` header rows
crash `create_report` with `ValueError: invalid literal for int()`.
A line starting with `#` is treated as a comment and is fine — and
matches CLAUDE.md's "BED-like outputs must have a `#`-prefixed
header" rule. Use:

```
#chrom	start	end	name
chr2	25246000	25259000	DNMT3A_2_promoter
```

Pair this with `--info-columns name` so the `name` column shows up in
the report's variant table.

## Track ordering

Render order is top-to-bottom in the viewer; put annotation FIRST so
gene tracks anchor the user's eye at the top, then per-sample BAM + 5mC
+ 5hmC triplets stacked below in sample-group order. The worked example
follows: gencode → EPDnew → CpGIslands → RepeatMasker → (per-sample:
BAM, 5mC, 5hmC).

## Colors (Okabe-Ito, group-paired)

For two-group studies (e.g., normal vs tumor) pick two color pairs out
of the Okabe-Ito palette so groups are pre-attentively distinguishable:

| Group  | 5mC color           | 5hmC color           |
|--------|---------------------|----------------------|
| Group A (normal) | `rgb(0,68,136)` blue | `rgb(204,121,167)` reddish-purple |
| Group B (tumor)  | `rgb(213,94,0)` vermillion | `rgb(230,159,0)` orange |

Annotation track colors (also Okabe-Ito): EPDnew = vermillion
`rgb(213,94,0)`, CpG islands = bluish-green `rgb(0,158,115)`,
RepeatMasker = sky-blue `rgb(86,180,233)`.

`scripts/generate_tracks_json.py` reads these from a `group_colors:`
map in the YAML spec, so a new group only needs one entry.

## EPDnew promoter track (hg38)

If your `databases_config.yaml` carries EPDnew for hg38, the suggested
keys are:

```yaml
reference_genomes:
  local:
    hg38:
      EPDnewCoding:    <path-to>/Hs_EPDnew.hg38.bed.gz
      EPDnewNonCoding: <path-to>/HsNC_EPDnew.hg38.bed.gz
```

Source: <https://epd.expasy.org/epd/human/human_database.php?db=human>

The skill driver doesn't load these by default — they're a methylation-
specific track. Either reference them directly from `tracks.json` or add
an `EPDnew` entry to a custom `annotation:` section in your YAML spec.
mm10 / mm39 / t2t builds don't ship with EPDnew.

## Reference-fasta vs `--genome hg38`

Always pass `--fasta` (skill driver default), never `--genome hg38`.
The igv.js bundled genome IDs require internet at view + render time;
`--fasta` + `--standalone` produces a fully-offline HTML. See
`references/best_practices.md` Reference section.

## When to use the apptainer SIF (mostly automatic)

The driver auto-detects whether to run via an apptainer SIF or PATH
`create_report` based on `SLURM_JOB_ID` and `$IGV_REPORTS_SIF`:

| Environment | Default | Why |
|---|---|---|
| Local / login node (`SLURM_JOB_ID` unset) | PATH `create_report` | No cold-start tax; simplest path. |
| Compute node under SLURM (`SLURM_JOB_ID` set) AND `$IGV_REPORTS_SIF` points at an existing SIF | apptainer + SIF | Fresh node = cold NFS cache = 1-2 M page faults on conda init (~2.5 us each). The SIF reads once into RAM, then stays warm. |
| Compute node BUT no SIF set/found | falls back to PATH `create_report` (logged) | Safe default; no surprise SIF-not-found error. |

Override either way with `--apptainer` / `--no-apptainer`. The decision
(auto vs. explicit) is logged at run start so post-mortems are unambiguous.

To set up the SIF once, pull from the Galaxy depot:

```bash
export IGV_REPORTS_SIF=/path/to/igv-reports_1.16.0.sif
wget -O "$IGV_REPORTS_SIF" \
  'https://depot.galaxyproject.org/singularity/igv-reports:1.16.0--pyh7cba7a3_0'
```

**Mandatory `--cleanenv` for the SIF (driver handles it).** Host RHEL 8
exports `SSL_CERT_FILE=/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem`
which doesn't exist inside the Galaxy-depot SIF. `create_report`'s
standalone build path makes an HTTPS GET (likely for the IGV.js
ideogram CDN) that crashes with `[SSL: CERTIFICATE_VERIFY_FAILED]`
mid-render. The driver always invokes `singularity exec --cleanenv ...`
to scrub host env vars before they enter the SIF, so users don't need
to remember the flag. If you call create_report from the SIF directly
(bypassing the driver), include `--cleanenv` yourself. See
`rules/apptainer_env_leak.md` for the full pattern.

## Worked example

`../examples/methylation_ont/` is the canonical end-to-end run:
- 4 COLO829 ONT samples (2 normal-blood × 2 tumor)
- 2 promoter windows (DNMT3A_2 + EZH2)
- 5mC + 5hmC bedGraph per sample (8 bedGraph files, pre-sliced)
- gencode + EPDnew + CGI + rmsk annotation slices

Run `bash examples/methylation_ont/build.sh` to regenerate the HTML;
read `examples/methylation_ont/recipe.md` for the slot-by-slot guide
to adapting it.

## Post-render verification

After building the HTML, run `scripts/verify_report.py` to confirm the
embedded content matches your inputs (region count, coordinates, track
names). For methylation viewers this catches the worst silent failure
mode — a render that succeeded for the wrong samples — which the input-
side validation alone can't catch.

```bash
python scripts/verify_report.py \
    --html         methylation_report.hg38.html \
    --sites        sites.hg38.bed \
    --track-config tracks.json \
    --min-size-mb  1.0 \
    --out          methylation_report.verify.tsv \
    --fail-on-fail
```

For `--track-config` builds the check uses the JSON's `name` fields; in
the YAML spec consumed by `generate_tracks_json.py`, those names are the
`name:` keys in `annotation:` and the auto-generated `<sample>`,
`<sample> 5mC`, `<sample> 5hmC` labels per sample. Picking specific
sample names in the YAML therefore drives the verifier's coverage —
generic names like "sample1" weaken the check.

**For cohort methylation runs** (multi-patient × per-sample HTMLs +
`index.html`), the cohort verifier (`scripts/verify_cohort.py`) is the
more relevant tool: it additionally catches sample-swap bugs (sample-2's
BAMs accidentally ending up in sample-1's HTML), missing samples, and
`index.html` drift. The methylation workflow is especially vulnerable to
sample-swap typos because each patient has multiple ONT runs with similar-
looking flowcell IDs (e.g., `PAU59807` vs `PAU61427`). Auto-invoked by
`build_igvreports.py --samplesheet`; see SKILL.md "Cohort-level
verification" for details.

## Cross-references

- `rules/igv.md` — bigwig-can't-be-sliced, y-axis-autoscale, UCSC
  `/dev/stdout` truncation; the rules that motivate this cheat-sheet.
- `rules/apptainer_vs_conda.md` — when the `--apptainer` flag pays off.
- `references/best_practices.md` — generic create_report flag reference;
  sites BED, tracks, reference, performance, pitfalls table.
- `examples/methylation_ont/recipe.md` — full slot-by-slot example doc.
- `CLAUDE.md` §3A — upstream ONT methylation pipeline (pod5 → dorado →
  modkit pileup → bedGraph + bigwig).
