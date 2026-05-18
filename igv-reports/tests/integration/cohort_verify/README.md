# cohort_verify_demo — regression test for `verify_cohort.py`

End-to-end check that the cohort verifier catches the four failure modes
it's designed to catch. Self-asserting — exit nonzero on any mismatch.

## What it does

1. Generates a 3-sample samplesheet (TSV) pointing at three real COLO829 ONT BAMs.
2. Calls `build_igvreports.py --samplesheet ... --no-verify` to produce
   `reports/sample_{1,2,3}.hg38.html` + `index.html`.
3. Runs `verify_cohort.py` against the clean cohort (expects all PASS).
4. Runs four corruption scenarios, each asserting the expected check FAILs:

   | Scenario | Corruption | Expected FAILs |
   |---|---|---|
   | A | Delete `sample_3.hg38.html` | `*/cohort_html_coverage`, `sample_3/html_exists` |
   | B | Replace `sample_1.hg38.html` with sample_2's content | `sample_1/sample_tracks_match`, `sample_1/no_cross_sample_contamination`, `sample_1/sample_id_embedded` |
   | C | Drop one `<li>` from `index.html` | `*/index_consistency` |
   | D | Truncate `sample_2.hg38.html` to 1 KB | `sample_2/html_min_size`, `sample_2/region_count` |

5. Cleans up generated `reports/`, samplesheet, sites BED, and logs/ on exit.

## Run

```bash
bash tests/integration/cohort_verify/scenarios.sh
```

Or as part of the full test suite:

```bash
bash tests/run_all.sh                  # all layers
bash tests/run_all.sh --integration-only
```

## BAM paths (parameterized)

Requires BAMs supplied via
env vars when running elsewhere:

```bash
IGV_REPORTS_TEST_BAM_1=/path/to/sample1.bam \
IGV_REPORTS_TEST_BAM_2=/path/to/sample2.bam \
IGV_REPORTS_TEST_BAM_3=/path/to/sample3.bam \
    bash tests/integration/cohort_verify/scenarios.sh
```

If a default doesn't exist and no env override is set, the script exits
**77** (POSIX skipped-test convention) and `run_all.sh` reports it as a
skip, not a failure.

Runtime: ~60-90 s on a warm node (3-sample cohort build at 1-bp point-variant
sites + 4 reverify cycles). Per-sample HTML ends up ~3-5 MB. Cold-cache
network reads of the underlying ONT BAMs can extend this to 2-3 min on
first invocation.

Disk: ~15 MB temporary under `reports/`, auto-cleaned via `trap`.

The sites BED uses 1-bp point-variant style coordinates (not 13 kb promoter
windows like the methylation example) so BAM slicing stays fast — we're
testing the verifier, not the renderer. Adapt for other workflows if you
want to exercise wider windows.

## How to provide BAMs

Set `IGV_REPORTS_TEST_BAM_{1,2,3}` to paths of three indexed BAMs you have
access to. The verifier doesn't care which BAMs — it only requires that the
three rows in the samplesheet declare *different* BAMs (so scenario B's
contamination check has signal). Without those env vars, the test exits 77
(POSIX skip).

## Why this is `integration`, not `smoke` or `unit`

This test depends on real BAMs and on `create_report` actually running, so
it can't fit in `tests/smoke/` (which uses only the committed COLO829 slice
fixture and runs in seconds) or `tests/unit/` (parser-only, no I/O).

For the parser-level regression checks that gave rise to this verifier,
see [tests/unit/test_verify_report.py](../../unit/test_verify_report.py).
