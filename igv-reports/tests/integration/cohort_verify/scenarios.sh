#!/usr/bin/env bash
# scenarios.sh — end-to-end integration test for scripts/verify_cohort.py.
#
# Builds a 3-sample cohort, snapshots a clean verify pass, then runs four
# corruption scenarios and asserts each triggers the expected check FAILs.
# Exit nonzero if any assertion misses.
#
# Runtime: ~6-8 min cold (cohort build dominates); ~30 s when cohort is cached.
# Disk: ~15 MB under ./reports/ (auto-cleaned on success unless KEEP_REPORTS=1).
#
# BAM source — three different indexed BAMs (any organism, any size). Defaults
# require env vars (no built-in defaults):
#   IGV_REPORTS_TEST_BAM_1, _2, _3
# Tests SKIP (exit 77) when defaults are unset and no override is provided.
set -euo pipefail

EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(cd "${EX_DIR}/../../.." && pwd)"
BUILD="${SKILL_DIR}/scripts/build_igvreports.py"
VERIFY="${SKILL_DIR}/scripts/verify_cohort.py"

# BAM sources — must be supplied via env vars. No defaults: integration tests
# require three indexed BAMs from a public ONT release (e.g. COLO829 from EBI).
BAM_S1="${IGV_REPORTS_TEST_BAM_1:-}"
BAM_S2="${IGV_REPORTS_TEST_BAM_2:-}"
BAM_S3="${IGV_REPORTS_TEST_BAM_3:-}"

for bam in "${BAM_S1}" "${BAM_S2}" "${BAM_S3}"; do
    if [[ -z "${bam}" || ! -f "${bam}" ]]; then
        echo "SKIP: integration test needs three indexed BAMs." >&2
        echo "      Set IGV_REPORTS_TEST_BAM_{1,2,3} to paths of three .bam files (each with sibling .bai)." >&2
        exit 77   # POSIX skipped-test convention
    fi
done

SHEET="${EX_DIR}/samplesheet.hg38.tsv"
SITES="${EX_DIR}/sites.hg38.bed"
OUTDIR="${EX_DIR}/reports"

cleanup() {
    if [[ -n "${KEEP_REPORTS:-}" ]]; then
        echo "(KEEP_REPORTS set — leaving artifacts in ${OUTDIR} for inspection)"
        return
    fi
    rm -rf "${OUTDIR}" "${SHEET}" "${SITES}" "${EX_DIR}/logs"
}
# Only cleanup on success — failures leave artifacts so they can be debugged.
trap 'rc=$?; if [[ $rc -eq 0 ]]; then cleanup; else echo "(scenarios.sh exited $rc — leaving artifacts in ${OUTDIR} for debug)"; fi' EXIT

# Requires `create_report` (pip install igv-reports) and `samtools` on PATH.
# If you use a conda env, activate it before invoking this script.

# --- 1. Generate fresh inputs --------------------------------------------------
# Point-variant style sites: 1-bp wide each, --flanking 300 = ~600 bp windows.
# Keeps BAM slicing fast (seconds) even with 100+ GB ONT BAMs and full
# annotation tracks. We're testing the verifier, not the renderer; tiny
# windows are sufficient.
cat >"${SITES}" <<EOF
#chrom	start	end	name
chr2	25246500	25246501	DNMT3A_SNV
chr7	148884000	148884001	EZH2_SNV
EOF

printf 'sample\tbam_tumor\tsites_bed\n' >"${SHEET}"
printf 'sample_1\t%s\t%s\n' "${BAM_S1}" "${SITES}" >>"${SHEET}"
printf 'sample_2\t%s\t%s\n' "${BAM_S2}" "${SITES}" >>"${SHEET}"
printf 'sample_3\t%s\t%s\n' "${BAM_S3}" "${SITES}" >>"${SHEET}"

# --- 2. Build cohort (3 HTMLs + index.html) -----------------------------------
# Skip rebuild if a complete cohort is already on disk (set REBUILD=1 to force).
# Lets you iterate on the verifier in seconds instead of waiting ~12 min to
# regenerate HTMLs that haven't changed.
if [[ -z "${REBUILD:-}" \
      && -f "${OUTDIR}/sample_1.hg38.html" \
      && -f "${OUTDIR}/sample_2.hg38.html" \
      && -f "${OUTDIR}/sample_3.hg38.html" \
      && -f "${OUTDIR}/index.html" ]]; then
    echo "=== reusing existing cohort in ${OUTDIR} (set REBUILD=1 to force rebuild) ==="
else
    echo "=== building cohort ==="
    python "${BUILD}" \
        --samplesheet "${SHEET}" \
        --genome hg38 \
        --flanking 300 \
        --type mutation \
        --info-columns name \
        --output-dir "${OUTDIR}" \
        --no-apptainer \
        --no-verify  # don't auto-verify during build — we exercise the verifier explicitly below
fi
echo

assert_status() {
    # assert_status <sample> <check> <expected_status> <verify_tsv>
    local sample="$1" check="$2" expected="$3" tsv="$4"
    local actual
    actual=$(awk -F'\t' -v s="$sample" -v c="$check" '$1==s && $2==c {print $3; exit}' "$tsv")
    if [[ "$actual" != "$expected" ]]; then
        echo "  FAIL ASSERTION: sample=$sample check=$check expected=$expected actual=${actual:-<missing>}"
        return 1
    fi
    echo "  OK   sample=$sample check=$check status=$actual"
}

# --- 3. Baseline verify (all PASS) --------------------------------------------
echo "=== scenario 0: baseline (all PASS expected) ==="
python "${VERIFY}" \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --out "${OUTDIR}/baseline.tsv" \
    --fail-on-fail >/dev/null
echo "  baseline: all PASS (verify exited 0)"
echo

# --- 4. Scenario A: missing HTML ----------------------------------------------
echo "=== scenario A: delete sample_3's HTML — C1 cohort_html_coverage should FAIL ==="
mv "${OUTDIR}/sample_3.hg38.html" "${OUTDIR}/sample_3.hg38.html.bak"
python "${VERIFY}" \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --out "${OUTDIR}/A.tsv" >/dev/null || true
assert_status "*"        "cohort_html_coverage" "FAIL" "${OUTDIR}/A.tsv"
assert_status "sample_3" "html_exists"          "FAIL" "${OUTDIR}/A.tsv"
mv "${OUTDIR}/sample_3.hg38.html.bak" "${OUTDIR}/sample_3.hg38.html"
echo

# --- 5. Scenario B: sample swap (sample_1.html now contains sample_2 data) ---
echo "=== scenario B: swap sample_1<-sample_2 — sample_tracks_match + id_embedded + contamination should FAIL on sample_1 ==="
cp "${OUTDIR}/sample_1.hg38.html" "${OUTDIR}/sample_1.hg38.html.bak"
cp "${OUTDIR}/sample_2.hg38.html" "${OUTDIR}/sample_1.hg38.html"
python "${VERIFY}" \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --out "${OUTDIR}/B.tsv" >/dev/null || true
assert_status "sample_1" "sample_tracks_match"          "FAIL" "${OUTDIR}/B.tsv"
assert_status "sample_1" "no_cross_sample_contamination" "FAIL" "${OUTDIR}/B.tsv"
assert_status "sample_1" "sample_id_embedded"           "FAIL" "${OUTDIR}/B.tsv"
mv "${OUTDIR}/sample_1.hg38.html.bak" "${OUTDIR}/sample_1.hg38.html"
echo

# --- 6. Scenario C: corrupt index.html ----------------------------------------
echo "=== scenario C: drop one <li> from index.html — C5 index_consistency should FAIL ==="
cp "${OUTDIR}/index.html" "${OUTDIR}/index.html.bak"
sed -i '/href="sample_2.hg38.html"/d' "${OUTDIR}/index.html"
python "${VERIFY}" \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --out "${OUTDIR}/C.tsv" >/dev/null || true
assert_status "*" "index_consistency" "FAIL" "${OUTDIR}/C.tsv"
mv "${OUTDIR}/index.html.bak" "${OUTDIR}/index.html"
echo

# --- 7. Scenario D: tiny HTML (truncation) ------------------------------------
echo "=== scenario D: truncate sample_2.html to 1 KB — html_min_size + parse failures expected ==="
cp "${OUTDIR}/sample_2.hg38.html" "${OUTDIR}/sample_2.hg38.html.bak"
head -c 1024 "${OUTDIR}/sample_2.hg38.html.bak" > "${OUTDIR}/sample_2.hg38.html"
python "${VERIFY}" \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --min-size-mb 1.0 \
    --out "${OUTDIR}/D.tsv" >/dev/null || true
assert_status "sample_2" "html_min_size"   "FAIL" "${OUTDIR}/D.tsv"
assert_status "sample_2" "region_count"    "FAIL" "${OUTDIR}/D.tsv"
mv "${OUTDIR}/sample_2.hg38.html.bak" "${OUTDIR}/sample_2.hg38.html"
echo

echo "=== all 4 scenarios PASSED — verify_cohort.py behaves as expected ==="
