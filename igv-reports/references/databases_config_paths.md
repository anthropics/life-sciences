# Databases-config YAML schema (for `--db-config` / `$IGV_REPORTS_DB_CONFIG`)

Optional. Without a databases YAML the driver still works — pass `--fasta` and
`--no-default-tracks` (plus any `--extra-track` you need) on every call.

The YAML is convenient when running across many regions/cohorts on the same
genome build: one file maps a short `--genome <id>` flag to the FASTA + the
default annotation tracks (CpG islands, gencode, RepeatMasker), so each
invocation stays short.

## Schema

```yaml
reference_genomes:
  local:
    <genome_id>:
      fasta:         <path>          # required
      gtf:           <path>          # gencode .gtf.gz or .gff3.gz (bgzip + tabix preferred)
      sizes:         <path>          # chrom.sizes (optional)
      CpGIslands:    <path>          # .bed (uncompressed or bgzip)
      repMaskerBed:  <path>          # .bed.gz (bgzip + tabix)
```

`<genome_id>` is the value you pass to `--genome`. Suggested IDs and aliases:

| `--genome` value     | YAML key            | Common alias |
|----------------------|---------------------|--------------|
| `hg38`               | `hg38`              | GRCh38       |
| `mm10`               | `mm10`              | GRCm38       |
| `mm39`               | `mm39`              | GRCm39       |
| `t2t` / `chm13`      | `t2t_CHM13v2_plusY` | T2T-CHM13v2  |
| `grch37` / `hg19`    | `GRCh37`            | hg19         |

The driver normalizes the input alias to the canonical YAML key. Extend
`GENOME_ALIASES` in `scripts/build_igvreports.py` if you need additional builds.

## Default-track resolution

For the `--genome` you pass, the driver tries to load three default tracks:

1. **CpG islands** → `CpGIslands` key
2. **Gene annotation** → `gtf` key (prefers a sibling `*.gff3.gz` if present)
3. **RepeatMasker** → `repMaskerBed` key

Any track absent from the YAML for that genome is logged as a warning and
skipped — the report still builds, just without that track.

## Gencode preference: GFF3 over GTF

If `gtf` points at `gencode.<version>.annotation.gtf.gz` and a sibling
`gencode.<version>.annotation.gff3.gz` exists in the same directory, the
driver prefers the GFF3 — it carries the full transcript / exon / CDS / UTR
detail that's most useful for read-level inspection at SV / fusion / integration
junctions. The GTF (gene-level) loads as a fallback.

Override with `--gencode-from-yaml` to force the YAML's `gtf` path regardless.

## EPDnew (methylation-specific)

`EPDnewCoding` / `EPDnewNonCoding` keys (BED.gz, bgzip + tabix) are
**not** auto-loaded — methylation-specific. Reference them explicitly via a
`--track-config tracks.json` entry when building a methylation viewer (see
`references/methylation_ont.md`).

## Missing tracks → workflow

1. Build or locate the BED / GFF3 / GTF.
2. If it needs bgzip + tabix conversion, run `scripts/prep_track.sh <path>`.
3. Add the path to your `databases_config.yaml` under the appropriate key, or
   pass it via `--extra-track <path>` for a one-off run.
