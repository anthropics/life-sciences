# tests/fixtures

Committed test fixtures derived from publicly released bioinformatics data.
Used by the smoke + integration test layers; safe to redistribute.

## tiny_colo829.hg38.bam

A 457 KB BAM (+ 85 KB `.bai`) sliced from Oxford Nanopore Technologies'
publicly released COLO829BL (matched normal) ONT reference dataset.

| Property | Value |
|---|---|
| **Source dataset** | ONT COLO829 / COLO829BL R10.4.1 5kHz sup basecalls |
| **ENA project** | PRJEB57425 |
| **Source flowcell** | PAU59807 (COLO829BL) |
| **Basecaller** | Dorado, model `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, `5mCG_5hmCG@latest,6mA@latest` |
| **Reference** | hg38 (`Homo_sapiens_assembly38.fasta`) |
| **Slice regions** | `chr2:25245000-25248000` (around DNMT3A), `chr7:148882000-148886000` (around EZH2) |
| **Subsample** | 20% reads, seed 42 (`samtools view --subsample 0.2 --subsample-seed 42`) |
| **Filtering** | `-F 1536` (drops PCR/optical dups + supplementary alignments — matches igv-reports' BamReader default) |
| **License** | The source data is openly released by ONT; this slice inherits that status. Slicing/subsampling is non-creative transformation. |

## Anchor sanity counts (used by smoke + integration tests)

| Region | `samtools view -c -F 1536` |
|---|---|
| `chr2:25246500-25246501` | **5** |
| `chr7:148884000-148884001` | **9** |

These counts are the contract: any change to the fixture (regeneration with
different params, etc.) must preserve these exact integers, or update the
constants in `tests/smoke/test_slice_count.py` and the integration `scenarios.sh`.

## Regenerate

```bash
bash tests/fixtures/build_fixtures.sh
```

Requires `samtools` (via PATH or `$SAMTOOLS_SIF`) and a local copy of the ONT
COLO829 release pointed to by `$COLO829BL_BAM`. Public source: ENA project
PRJEB57425.

## Why these regions

The two sites are coding mutations in well-known cancer driver genes
(DNMT3A R882, EZH2 Y646) at coordinates the demos already use. Picking
real loci keeps the test data biologically interpretable and lets the same
fixture exercise both the parser layer (anchors named for real variants
read naturally) and the slice-decode layer (read counts you can sanity-check
in IGV against the source BAM if needed).
