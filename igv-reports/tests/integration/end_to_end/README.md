# tests/integration/end_to_end

End-to-end smoke test against the **committed** `tests/fixtures/tiny_colo829.hg38.bam`
fixture (457 KB). Unlike the other integration scenarios (`anchor_verify`,
`cohort_verify`), this one needs no shared-storage access and runs in
~30 s — so it ships in CI.

## What it exercises (not via mocks)

1. **`build_igvreports.py --bam ...`** actually invokes `create_report`
   against the fixture BAM with a synthesized minimal FASTA. Produces a
   real ~2 MB HTML.
2. **`verify_report.py`** parses the HTML's `tableJson` + `sessionDictionary`,
   confirms region count + track presence.
3. **`verify_anchors.py generate`** counts reads in the source BAM at the
   three sites; asserts the counts match the frozen contract documented in
   `tests/fixtures/README.md` (`chr2:25246500-25246501 = 5`,
   `chr7:148884000-148884001 = 9`).
4. **`verify_anchors.py verify`** decodes the embedded BAM slices from the
   freshly-built HTML and confirms the same counts — closes the loop on the
   create_report ↔ source-BAM round trip.
5. **`--also-png` (optional)** runs the same pipeline with the PNG sidecar
   path. SKIPs cleanly when `igver` isn't installed or fails (the
   documented silent-failure mode in `rules/igv.md`).

## What it catches that unit tests don't

- `create_report` flag rename / removal on upstream version bumps
- HTML structural changes from upstream (e.g. session-dict layout drift)
- Driver regressions on the **non-mock** code path
- Off-MSKCC portability bugs — the test runs against the committed fixture
  with no `shared-storage` dependency, so CI exercises the same code paths
  external users would hit

## Runtime

| Step | Duration |
|---|---|
| `create_report` (3 regions × 1 BAM × 300 bp flanking) | ~2 s |
| structural `verify_report.py` | <1 s |
| `verify_anchors.py generate` + `verify` | ~5 s |
| `--also-png` (if igver available) | ~5 s |
| **total** | **~14 s** |

## Prereqs

- `create_report` on PATH (`pip install -U 'igv-reports>=1.16.0'`)
- `samtools` on PATH (provided by the smoke layer prereqs)
- `python3` on PATH

If `create_report` is missing the test exits 77 (skipped) rather than
failing — same convention as the other integration scenarios.

## Knobs

- `KEEP_REPORTS=1` — leave the `out/` directory in place after a successful
  run for manual inspection.
- `IGV_REPORTS_PY=/path/to/python` — pin the python interpreter (the
  default search is conda's snakemake env → `python3` on PATH).
- `IGVER_CMD='apptainer exec /path/to/igver.sif igver'` — provide a working
  `igver` invocation so step 8 (`--also-png`) actually exercises the PNG
  pipeline rather than SKIPping.
