# igv-reports best practices

Authoritative companion to the skill. Read this when something fails in a
way the SKILL.md pitfalls table doesn't cover, or when introducing a new
input format / track type.

## Sites/regions input

Supported by `create_report`:
- **VCF** — variant table is built from CHROM/POS/ID/REF/ALT plus any
  `--info-columns` you surface from INFO and `--sample-columns` from
  FORMAT. Use `--idlink "https://url/$$"` to make ID a clickable link.
- **BED** — fields parsed by position: `chr / start / end [/ name]`.
  A **non-comment header row** (e.g., `chrom start end name`) crashes
  `create_report` with `ValueError: invalid literal for int()` because
  the parser tries to `int()` the string `start`. A `#`-prefixed comment
  header (e.g., `#chrom\tstart\tend\tname`) IS accepted — `create_report`
  skips lines starting with `#`. This matches the lab's "BED-like outputs
  must have a `#`-prefixed header" convention in CLAUDE.md.
- **MAF** — Mutation Annotation Format (TCGA standard).
- **BEDPE** — paired-end / fusion / SV format. With `--type fusion` each
  row is rendered as a multi-locus split-screen view.
- **Generic TSV** — any tab-delimited file. Requires `--sequence`,
  `--begin`, `--end` to name the chrom/start/end columns. Add
  `--zero_based` if 0-based.

**File-extension dispatch**: igv-reports picks the parser by extension,
not content. `.bed` → BED parser (which IGNORES `--sequence/--begin/--end`).
If you want a TSV-with-header parsed by name, the extension must NOT be
`.bed`/`.vcf`/`.gff3`/`.maf` — use `.tsv` or `.txt`.

The project's `enforce-genome-tag.sh` hook requires a genome tag in the
filename: `sites.hg38.bed`, not `sites.bed`.

## Tracks

Supported track formats: BAM, CRAM, VCF, BED, GFF3, GTF, WIG, BEDGRAPH.

**Indexing**:
- BAM/CRAM/VCF MUST be indexed (`.bai`/`.crai`/`.tbi` sidecar).
- Large `.bed.gz` / `.gff3.gz` / `.gtf.gz` SHOULD be tabix-indexed
  (`.tbi` sidecar) and **must be true bgzip** — not plain gzip.
- Check format with `file <name>` — true bgzip says
  `gzip compressed data, extra field, original size 0`. Plain gzip
  has no "extra field". igv-reports trips on plain-gzip .gff3.gz with
  cryptic `UnicodeDecodeError: byte 0x8b at position 1` — that 0x8b is
  the gzip magic byte the parser is reading as text.

**Sortedness**: gencode and many other GFF/GTF distributions interleave
records by feature type at the same locus (gene → transcript → exon → CDS →
exon → CDS → ...) rather than strictly position-sorted within each
chromosome. tabix requires pos-sorted within chr. Fix:
`sort -k1,1 -k4,4n` on the body, then bgzip + tabix. The `prep-track`
script in this skill does the full pipeline with backup.

**Track render order**: the order you pass to `--tracks` is the order
they appear in the IGV.js viewer (top-to-bottom). Convention:
1. BAM/CRAM (the data you want to evaluate)
2. VCF (the calls being inspected)
3. Annotation tracks (genes, regulatory, repeats, CGI)

The skill defaults always render annotation tracks LAST so they sit at
the bottom and don't push the read evidence off-screen.

## Reference

One of `--fasta`, `--twobit`, or `--genome` is required.

- `--fasta /path/to/local.fa` (with `.fai`) — fully offline, supports
  custom or combined references (e.g., host + viral).
- `--genome hg38` — uses igv.js bundled IDs, but **requires internet at
  view AND render time** because igv.js fetches the bundled genome.
  Avoid for HPC/offline.
- `--twobit` — alternative reference in 2bit format.

For combined viral+host references, the single FASTA must include all
contigs, and any per-contig tracks must align (e.g., HTLV1_features.bed
must use the same contig name as in the FASTA).

## Window sizing

`--flanking N` (igv-reports default 1000, this skill default **300**)
adds N bp on either side of each site.

| Use case | Recommended flanking |
|---|---|
| Point variants (SNV/indel) | 50–200 bp |
| SV / integration breakpoints | 300–1000 bp (this skill: 300) |
| Whole-gene context | gene length + 5–10 kb |

`--maxlen N` (default 10,000) — variants exceeding this length switch to
split-screen multilocus view automatically. Useful for SVs > 10 kb.

`--window N` — initial visible window inside the embedded igv.js viewer
(if not supplied, igv.js defaults to 41 bp, which is too narrow for
read-level inspection). Set to ~`2 × flanking` so the user lands on the
full embedded slice.

## Output

- `--standalone` embeds all igv.js JS in the HTML → fully offline,
  4–11 MB per patient typical for cohort runs.
- `--no-embed` keeps external URLs → smaller HTML but online required.
  Avoid for HPC/sharing-by-email.

Per-region BAM data is ALWAYS sliced and embedded by default; only the
flanking-sized portion of large BAMs ships in the HTML — so the HTML stays
manageable even when input BAMs are 100+ GB.

## Variant table customization

For VCF input:
- `--info-columns SVTYPE SVLEN ALIGNED_POS DR DV VAF` surfaces those
  INFO fields as table columns.
- `--info-columns-prefixes ANN_ HTLV1_` includes any INFO field starting
  with the listed prefixes.
- `--sample-columns DP AD GT` (with optional `--samples NAME`) surfaces
  per-sample FORMAT fields.
- `--idlink "https://example.com/$$"` makes the VCF ID column clickable
  with `$$` replaced by the ID value.

Order of operations: include `--info-columns` for the call-quality fields
your reviewer needs to see at a glance; the rest is one click into the
variant detail.

## Performance / size control

- `--subsample 0.0-1.0` — keep a fraction of BAM alignments per region.
  Use for very deep BAMs (>100×) where the rendered viewer would be
  read-cluttered.
- `--exclude-flags 1536` (default) — excludes duplicates and QC-fail
  reads. Set to 0 to keep everything.
- Render time scales roughly linearly with `n_regions × n_tracks`. The
  ATLL cohort run (6 patients × 1–3 integrations + HTLV1 + EBV regions,
  6 tracks) took ~2 min/patient with the gene-level GFF and ~3 min/patient
  with the full annotation.

## Pitfalls observed in production

| Symptom | Root cause | Fix |
|---|---|---|
| `ValueError: invalid literal for int() with base 10: 'start'` | Non-comment header row in BED sites file | Prefix the header with `#` (skipped by create_report and matches lab convention); or strip it entirely |
| `UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b` | igv-reports reading bgzip as text (file actually plain-gzip but with `.gz` ext) | Convert with prep-track; verify with `file <name>` |
| `tabix: not BGZF` | Plain gzip masquerading as `.gz` | `gunzip → bgzip → tabix` |
| `tabix: out of order` | GFF/GTF/BED records not pos-sorted within chr | `sort -k1,1 -k4,4n` first |
| Empty annotation track in viewer | Tabix lookup returns nothing in window. Often correct biology (e.g., CGI-distal site) — verify with `tabix file region` |
| Title shows weird characters | Unicode em-dash (`—`) in `--title` got mangled by shell escaping | Use plain ASCII `-` |
| HTML loads but viewer is blank | `--genome hg38` without internet at view time | Use `--fasta` + `--standalone` |
| `tabix` index missing for a track | igv-reports looked for `<track>.tbi`, not present | Re-run `tabix -p <gff|gtf|bed>` |
| `samtools index` errors mid-render | BAM index stale (BAM modified after `.bai`) | `samtools index -@ 4 file.bam` |
| Output HTML size much larger than expected | Some region accidentally spans Mb-scale (e.g., a row with start=0 end=chrom_length); flanking compounds this | Validate the sites BED — `awk '$3-$2 > 1e6'` to find offenders |

## See also

- Official docs: https://github.com/igvteam/igv-reports
- igv.js track config schema: https://github.com/igvteam/igv.js/wiki/Tracks-2.0
- This skill's `references/databases_config_paths.md` for which YAML keys
  hold which tracks per genome.
