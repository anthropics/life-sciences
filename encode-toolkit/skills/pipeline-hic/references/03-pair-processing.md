# Pair Processing with pairtools

pairtools is the ENCODE-standard tool for Hi-C pair classification,
filtering, and deduplication. It processes aligned BAMs into sorted,
deduplicated .pairs files ready for matrix generation.

## Parse: Classify Read Pairs

Convert aligned BAM to .pairs format with pair type classification:

```bash
# pairtools sort hands --tmpdir to GNU sort, which does not create it:
# the directory must already exist.
mkdir -p tmp

pairtools parse \
    --chroms-path chrom.sizes \
    --min-mapq 30 \
    --walks-policy mask \
    --max-inter-align-gap 30 \
    --nproc-in 4 \
    --nproc-out 4 \
    --output-stats sample.parse_stats.txt \
    sample.paired.bam \
    | pairtools sort \
        --nproc 4 \
        --tmpdir $PWD/tmp \
        -o sample_parsed_sorted.pairs.gz
```

`main.nf` publishes the statistics file as `pairs/{sample}.parse_stats.txt`.
It holds the pair-type breakdown for the whole library, which is the only
place those counts are available (see below).

### Parse Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--min-mapq 30` | MAPQ 30 | Filter low-confidence alignments |
| `--walks-policy mask` | Mask walks | Handle complex ligation events |
| `--max-inter-align-gap 30` | 30 bp | Maximum gap between split alignments |

### Pair Types Output

pairtools assigns each pair a two-letter code, one letter per side: U unique,
R rescued, M multi, N null (unmapped), W walk, D duplicate, X corrupt.

| Code | Meaning | Use |
|------|---------|-----|
| UU | Both uniquely mapped | Primary contacts |
| UR/RU | One unique, one rescued | Valid with caution |
| MU | One multi-mapped, one unique | Ambiguous, excluded |
| MM | Both multi-mapped | Excluded |
| NU | One unique, one unmapped | Excluded |
| NM | One unmapped, one multi-mapped | Excluded |
| NN | Both unmapped | Excluded |
| WW | Complex walk (multiple ligations), masked by `--walks-policy mask` | Excluded |
| DD | Duplicate | Removed in dedup step |
| XX | Corrupt record | Excluded |

## Sort Pairs

Pairs must be sorted by genomic position for deduplication:

```bash
pairtools sort \
    --nproc 4 \
    --tmpdir /tmp/ \
    sample_parsed.pairs.gz \
    -o sample_sorted.pairs.gz
```

## Deduplicate

Remove PCR/optical duplicates based on alignment positions:

```bash
pairtools dedup \
    --nproc-in 4 \
    --nproc-out 4 \
    --mark-dups \
    --output-stats sample.dedup_stats.txt \
    -o sample.dedup.pairs.gz \
    sample_sorted.pairs.gz
```

Both files are published: `pairs/{sample}.dedup.pairs.gz` and
`pairs/{sample}.dedup_stats.txt`.

### Dedup Statistics

The stats file reports:
- Total pairs processed
- Unique pairs retained
- PCR duplicate pairs removed
- Optical duplicate pairs removed
- Complexity estimate

Expected duplication rate: 10-40% depending on library complexity and depth.

## Filter for Valid Contacts

Select only UU pairs for contact matrix generation:

```bash
pairtools select \
    '(pair_type == "UU")' \
    sample.dedup.pairs.gz \
    -o sample_valid.pairs.gz
```

This is what `main.nf` does; the selected pairs file itself is an intermediate
and is not published.

For higher sensitivity (at cost of some noise), you can include rescued pairs
manually. The workflow does not offer this as an option:

```bash
pairtools select \
    '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
    sample.dedup.pairs.gz \
    -o sample_valid_rescued.pairs.gz
```

## Pair Statistics

Generate detailed contact statistics:

```bash
pairtools stats \
    sample_valid.pairs.gz \
    -o sample_contact_stats.txt
```

`main.nf` runs exactly this on the UU-selected pairs and publishes it as
`qc/{sample}.contact_stats.txt`.

Key metrics from the stats output:
- **cis contacts**: Same chromosome
- **trans contacts**: Different chromosomes
- **cis >20kb**: Long-range cis contacts (biologically meaningful)
- **cis <20kb**: Short-range, often ligation artifacts

Because the input is UU-only, the pair-type distribution in this file is 100%
UU by construction. Read pair types from `pairs/{sample}.parse_stats.txt`
instead.

## Cis/Trans Ratio

The cis/trans ratio is a key QC metric. Match whole keys: `pairtools stats`
also emits `cis_1kb+` ... `cis_40kb+` rows, and a regex like `/cis/` would
pick up the last of those instead of the `cis` total.

```bash
awk '$1=="cis" {c=$2} $1=="trans" {t=$2} END {
    print "Cis:", c;
    print "Trans:", t;
    print "Cis/Trans ratio:", c/t
}' sample_contact_stats.txt
```

For thresholds, use the single QC table in `SKILL.md` (cis/trans >1.5 pass,
1.0-1.5 warning, <1.0 fail).
