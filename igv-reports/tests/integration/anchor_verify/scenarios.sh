#!/usr/bin/env bash
# scenarios.sh — end-to-end integration test for scripts/verify_anchors.py.
#
# Builds a 2-sample cohort, freezes BAM-read-count anchors from the source
# BAMs, verifies the clean cohort, then runs four corruption scenarios and
# asserts each triggers the expected PASS / FAIL / SKIP outcomes.
#
# Runtime: ~6-8 min cold (cohort build dominates); ~15 s when cohort is cached.
# Disk: ~10 MB under ./reports/ (auto-cleaned on success unless KEEP_REPORTS=1).
#
# BAM source — two different indexed BAMs (any organism, any size). Defaults
# require env vars (no built-in defaults):
#   IGV_REPORTS_TEST_BAM_1, _2
# Tests SKIP (exit 77) when defaults are unset and no override is provided.
set -euo pipefail

EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(cd "${EX_DIR}/../../.." && pwd)"
BUILD="${SKILL_DIR}/scripts/build_igvreports.py"
ANCHORS="${SKILL_DIR}/scripts/verify_anchors.py"

# BAM sources — must be supplied via env vars. No defaults: integration tests
# require two indexed BAMs from a public release.
BAM_S1="${IGV_REPORTS_TEST_BAM_1:-}"
BAM_S2="${IGV_REPORTS_TEST_BAM_2:-}"

for bam in "${BAM_S1}" "${BAM_S2}"; do
    if [[ -z "${bam}" || ! -f "${bam}" ]]; then
        echo "SKIP: integration test needs two indexed BAMs." >&2
        echo "      Set IGV_REPORTS_TEST_BAM_{1,2} to paths of two .bam files (each with sibling .bai)." >&2
        exit 77   # POSIX skipped-test convention
    fi
done

SHEET="${EX_DIR}/samplesheet.hg38.tsv"
SITES="${EX_DIR}/sites.hg38.bed"
OUTDIR="${EX_DIR}/reports"
ANCHORS_TSV="${EX_DIR}/anchors.hg38.tsv"

cleanup() {
    if [[ -n "${KEEP_REPORTS:-}" ]]; then
        echo "(KEEP_REPORTS set — leaving artifacts in ${OUTDIR} and ${EX_DIR}/anchors* for inspection)"
        return
    fi
    rm -rf "${OUTDIR}" "${SHEET}" "${SITES}" "${ANCHORS_TSV}" \
           "${EX_DIR}/anchors.corrupted.tsv" "${EX_DIR}/anchors.min.tsv" \
           "${EX_DIR}/anchors.subset.tsv" "${EX_DIR}/logs"
}
trap 'rc=$?; if [[ $rc -eq 0 ]]; then cleanup; else echo "(scenarios.sh exited $rc — leaving artifacts for debug)"; fi' EXIT

# Requires `create_report` (pip install igv-reports) and `samtools` on PATH.
# If you use a conda env, activate it before invoking this script.

# Defensive restore: a previous run may have died mid-corruption leaving
# .bak files. Bring HTMLs back to their original state before we start.
for f in "${OUTDIR}"/*.hg38.html.bak; do
    [[ -f "$f" ]] && mv "$f" "${f%.bak}" && echo "(restored ${f%.bak} from .bak)"
done 2>/dev/null || true

# --- 1. Inputs -----------------------------------------------------------------
# Two SNV-style point sites; --flanking 300 keeps BAM slicing in seconds even
# at 167 GB source BAMs. We're testing the verifier, not the renderer.
cat >"${SITES}" <<EOF
#chrom	start	end	name
chr2	25246500	25246501	DNMT3A_SNV
chr7	148884000	148884001	EZH2_SNV
EOF

printf 'sample\tbam_tumor\tsites_bed\n'  >"${SHEET}"
printf 'sample_1\t%s\t%s\n' "${BAM_S1}" "${SITES}" >>"${SHEET}"
printf 'sample_2\t%s\t%s\n' "${BAM_S2}" "${SITES}" >>"${SHEET}"

# --- 2. Generate anchors from source BAMs --------------------------------------
echo "=== generate: freezing samtools-view counts as anchors ==="
python "${ANCHORS}" generate \
    --samplesheet "${SHEET}" \
    --sites "${SITES}" \
    --out "${ANCHORS_TSV}" 2>&1 | tail -6
echo

# --- 3. Build cohort -----------------------------------------------------------
if [[ -z "${REBUILD:-}" \
      && -f "${OUTDIR}/sample_1.hg38.html" \
      && -f "${OUTDIR}/sample_2.hg38.html" ]]; then
    echo "=== reusing existing cohort in ${OUTDIR} (set REBUILD=1 to force) ==="
else
    echo "=== building 2-sample cohort (this takes ~5-7 min on warm node) ==="
    python "${BUILD}" \
        --samplesheet "${SHEET}" \
        --genome hg38 \
        --flanking 300 \
        --type mutation \
        --info-columns name \
        --output-dir "${OUTDIR}" \
        --no-apptainer \
        --no-verify  # auto-verify is structural; we exercise the anchor verifier ourselves below
fi
echo

assert_status() {
    # assert_status <sample> <region> <expected_status> <verify_tsv>
    local sample="$1" region="$2" expected="$3" tsv="$4"
    local actual
    actual=$(awk -F'\t' -v s="$sample" -v r="$region" '$1==s && $3==r {print $4; exit}' "$tsv")
    if [[ "$actual" != "$expected" ]]; then
        echo "  FAIL ASSERTION: sample=$sample region=$region expected=$expected actual=${actual:-<missing>}"
        return 1
    fi
    echo "  OK   sample=$sample region=$region status=$actual"
}

# --- 4. Scenario 0: clean cohort, all PASS -------------------------------------
echo "=== scenario 0: clean — all anchors expected PASS ==="
python "${ANCHORS}" verify-cohort \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --anchors "${ANCHORS_TSV}" \
    --out "${OUTDIR}/scenario0.tsv" \
    --fail-on-fail >/dev/null
echo "  baseline: 4/4 PASS (verify-cohort exited 0)"
echo

# --- 5. Scenario A: tolerance violation ----------------------------------------
echo "=== scenario A: corrupt expected count outside tolerance — FAIL on diff_ratio ==="
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next} NR==2 {$6=9999; print; next} {print}' "${ANCHORS_TSV}" > "${EX_DIR}/anchors.corrupted.tsv"
python "${ANCHORS}" verify-cohort \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --anchors "${EX_DIR}/anchors.corrupted.tsv" \
    --out "${OUTDIR}/A.tsv" >/dev/null || true
assert_status "sample_1" "chr2:25246500-25246501" "FAIL" "${OUTDIR}/A.tsv"
assert_status "sample_1" "chr7:148884000-148884001" "PASS" "${OUTDIR}/A.tsv"
echo

# --- 6. Scenario B: min/max bound violation ------------------------------------
echo "=== scenario B: anchor min=1000 (real count ~56) — FAIL on min ==="
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next} NR==2 {$8=1000; print; next} {print}' "${ANCHORS_TSV}" > "${EX_DIR}/anchors.min.tsv"
python "${ANCHORS}" verify-cohort \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --anchors "${EX_DIR}/anchors.min.tsv" \
    --out "${OUTDIR}/B.tsv" >/dev/null || true
assert_status "sample_1" "chr2:25246500-25246501" "FAIL" "${OUTDIR}/B.tsv"
echo

# --- 7. Scenario C: corrupt data URL inside HTML — FAIL on decode --------------
echo "=== scenario C: mangle a session's base64 payload — FAIL on session decode ==="
cp "${OUTDIR}/sample_1.hg38.html" "${OUTDIR}/sample_1.hg38.html.bak"
# Replace one base64 chunk inside a session data URL. The H4sI prefix is the
# base64-encoded gzip magic 0x1f 0x8b 0x08; mangling it breaks the gunzip step
# that decodes the session, simulating arbitrary HTML tampering.
sed -i 's|data:application/gzip;base64,H4sI|data:application/gzip;base64,XXXX|g' "${OUTDIR}/sample_1.hg38.html"
python "${ANCHORS}" verify-cohort \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --anchors "${ANCHORS_TSV}" \
    --out "${OUTDIR}/C.tsv" >/dev/null || true
# Both regions in sample_1 should FAIL (sed hits every session URL in the file).
assert_status "sample_1" "chr2:25246500-25246501" "FAIL" "${OUTDIR}/C.tsv"
assert_status "sample_1" "chr7:148884000-148884001" "FAIL" "${OUTDIR}/C.tsv"
# sample_2 unaffected.
assert_status "sample_2" "chr2:25246500-25246501" "PASS" "${OUTDIR}/C.tsv"
mv "${OUTDIR}/sample_1.hg38.html.bak" "${OUTDIR}/sample_1.hg38.html"
echo

# --- 8. Scenario D: anchor missing for a (sample, region) — SKIP not FAIL ------
echo "=== scenario D: drop sample_1's chr2 anchor — that region SKIPs, others PASS ==="
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next} !($1=="sample_1" && $3=="chr2"){print}' "${ANCHORS_TSV}" > "${EX_DIR}/anchors.subset.tsv"
python "${ANCHORS}" verify-cohort \
    --samplesheet "${SHEET}" \
    --reports-dir "${OUTDIR}" \
    --genome hg38 \
    --anchors "${EX_DIR}/anchors.subset.tsv" \
    --out "${OUTDIR}/D.tsv" \
    --fail-on-fail >/dev/null
# The dropped anchor shouldn't appear at all (nothing to verify). Remaining anchors PASS.
n_rows=$(awk -F'\t' 'NR>1 && $1=="sample_1" && $3=="chr2"' "${OUTDIR}/D.tsv" | wc -l)
if [[ "${n_rows}" -ne 0 ]]; then
    echo "  FAIL ASSERTION: sample_1/chr2 should NOT appear (dropped anchor) but got ${n_rows} rows"
    exit 1
fi
echo "  OK   sample_1/chr2 anchor dropped — no row emitted"
assert_status "sample_1" "chr7:148884000-148884001" "PASS" "${OUTDIR}/D.tsv"
assert_status "sample_2" "chr2:25246500-25246501" "PASS" "${OUTDIR}/D.tsv"
echo

echo "=== all 4 scenarios PASSED — verify_anchors.py behaves as expected ==="
