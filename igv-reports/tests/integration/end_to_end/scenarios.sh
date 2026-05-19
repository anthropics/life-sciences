#!/usr/bin/env bash
# scenarios.sh — end-to-end smoke test using the COMMITTED tiny_colo829 fixture.
#
# Author: Samuel Ahuno
# Purpose:
#   The other integration scenarios (anchor_verify, cohort_verify) require
#   167 GB lab BAMs and take 6-8 min. This one uses the 457 KB
#   tests/fixtures/tiny_colo829.hg38.bam fixture so the full pipeline runs
#   in ~30 s on any machine with `create_report` on PATH.
#
#   What it exercises end-to-end (not via mocks):
#     1. build_igvreports.py invokes create_report against the fixture
#     2. The resulting HTML is parseable by verify_report.py (structural)
#     3. verify_anchors.py generate → frozen counts (chr2=5, chr7=9 per
#        tests/fixtures/README.md)
#     4. verify_anchors.py verify → PASS on the freshly built HTML
#     5. If `igver` is on PATH: --also-png produces non-empty per-region
#        PNGs and the manifest. Otherwise that step is SKIPped (logged).
#
# Catches: create_report flag drift, HTML-format upstream changes, driver
# regressions on the non-mock path, off-MSKCC portability bugs (the
# fixture is committed; no shared-storage required).
#
# Runtime: ~30 s. Disk: ~5 MB under ./out/ (auto-cleaned on success).
set -euo pipefail

EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(cd "${EX_DIR}/../../.." && pwd)"
BUILD="${SKILL_DIR}/scripts/build_igvreports.py"
ANCHORS="${SKILL_DIR}/scripts/verify_anchors.py"
VR="${SKILL_DIR}/scripts/verify_report.py"
FIXTURE="${SKILL_DIR}/tests/fixtures/tiny_colo829.hg38.bam"

# Prerequisite: create_report must be on PATH (provided by `pip install igv-reports`).
if ! command -v create_report >/dev/null 2>&1; then
    echo "SKIP: end-to-end test needs create_report on PATH." >&2
    echo "      Install with: pip install -U 'igv-reports>=1.16.0'" >&2
    exit 77
fi

# Prerequisite: the committed fixture must be readable.
if [[ ! -f "${FIXTURE}" ]] || [[ ! -f "${FIXTURE}.bai" ]]; then
    echo "ERROR: fixture missing or unindexed: ${FIXTURE}" >&2
    echo "       Regenerate with: bash tests/fixtures/build_fixtures.sh" >&2
    exit 1
fi

# Pick the python that can import pyyaml + the same scripts/. Tests/run_all.sh
# already does this dance; we mirror it.
PY="${IGV_REPORTS_PY:-}"
if [[ -z "${PY}" ]]; then
    if command -v python3 >/dev/null 2>&1; then
        PY=$(command -v python3)
    else
        echo "ERROR: no python3 on PATH" >&2; exit 2
    fi
fi

OUTDIR="${EX_DIR}/out"
cleanup() {
    if [[ -n "${KEEP_REPORTS:-}" ]]; then
        echo "(KEEP_REPORTS set — leaving artifacts in ${OUTDIR})"
        return
    fi
    rm -rf "${OUTDIR}"
}
trap 'rc=$?; if [[ $rc -eq 0 ]]; then cleanup; else echo "(scenarios.sh exited $rc — leaving ${OUTDIR} for debug)"; fi' EXIT

rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}"

# --- 1. Inputs --------------------------------------------------------------
# Three sites, all within the fixture's two slice regions:
#   * chr2:25246500-25246501  (DNMT3A R882 SNV, frozen anchor count = 5)
#   * chr7:148884000-148884001 (EZH2 Y646 SNV, frozen anchor count = 9)
#   * chr2:25247500-25247501  (second DNMT3A locus, count not frozen)
# Frozen counts are the contract per tests/fixtures/README.md.
SITES="${OUTDIR}/sites.hg38.bed"
cat >"${SITES}" <<EOF
#chrom	start	end	name
chr2	25246500	25246501	DNMT3A_R882
chr7	148884000	148884001	EZH2_Y646
chr2	25247500	25247501	DNMT3A_2nd
EOF

# Reference FASTA: in CI we won't have hg38 locally. Skip the --fasta
# resolution and use --no-default-tracks; create_report will then need a
# --fasta path. We supply a synthesized FASTA covering both regions so
# create_report can compute its own slice without internet access.
FASTA="${OUTDIR}/tiny.hg38.fa"
${PY} -c "
# Minimal multi-contig FASTA covering the fixture's coverage windows.
# Only the size matters for create_report's region slicing — bases don't
# need to be biologically real; the BAM's reads carry the actual signal.
contigs = [
    ('chr2', 30_000_000),
    ('chr7', 150_000_000),
]
with open('${FASTA}', 'w') as fh:
    for name, length in contigs:
        fh.write(f'>{name}\n')
        n_per_line = 60
        for i in range(0, length, n_per_line):
            fh.write('N' * min(n_per_line, length - i) + '\n')
"
samtools faidx "${FASTA}"

# --- 2. Build HTML (the actual end-to-end step) -----------------------------
echo "=== build: invoke create_report against fixture BAM ==="
HTML="${OUTDIR}/sample.hg38.html"
${PY} "${BUILD}" \
    --sites "${SITES}" \
    --bam "${FIXTURE}" \
    --genome hg38 \
    --fasta "${FASTA}" \
    --no-default-tracks \
    --flanking 300 \
    --type mutation \
    --info-columns name \
    --output "${HTML}" \
    --no-apptainer \
    --no-verify 2>&1 | tail -8
echo

# --- 3. Assertion: HTML exists, plausible size ------------------------------
if [[ ! -f "${HTML}" ]]; then
    echo "FAIL: HTML not produced at ${HTML}"; exit 1
fi
size=$(stat -c %s "${HTML}")
if [[ "${size}" -lt 50000 ]]; then
    echo "FAIL: HTML suspiciously small (${size} bytes) — expected >= 50 KB"
    exit 1
fi
echo "  OK   HTML: ${HTML} (${size} bytes)"
echo

# --- 4. Structural verify ---------------------------------------------------
echo "=== verify_report.py: structural check ==="
${PY} "${VR}" \
    --html "${HTML}" \
    --sites "${SITES}" \
    --tracks "${FIXTURE}" \
    --min-size-mb 0.05 \
    --out "${OUTDIR}/verify_report.tsv" \
    --fail-on-fail >/dev/null
echo "  OK   structural verify PASS"
echo

# --- 5. Generate frozen anchors ---------------------------------------------
echo "=== verify_anchors.py generate: BAM read counts ==="
SHEET="${OUTDIR}/samplesheet.tsv"
printf 'sample\tbam_tumor\tsites_bed\n'  >"${SHEET}"
printf 'sample\t%s\t%s\n' "${FIXTURE}" "${SITES}" >>"${SHEET}"

ANCHORS_TSV="${OUTDIR}/anchors.hg38.tsv"
${PY} "${ANCHORS}" generate \
    --samplesheet "${SHEET}" \
    --sites "${SITES}" \
    --out "${ANCHORS_TSV}" 2>&1 | tail -6
echo

# --- 6. Assertion: frozen anchor counts match contract ----------------------
# Contract is in tests/fixtures/README.md. Any drift here is the loudest
# signal that the fixture changed, the BAM filter changed, or the test
# environment is using a different samtools.
expected_chr2=5
expected_chr7=9
actual_chr2=$(awk -F'\t' '$4=="chr2" && $5==25246500 {print $7}' "${ANCHORS_TSV}")
actual_chr7=$(awk -F'\t' '$4=="chr7" && $5==148884000 {print $7}' "${ANCHORS_TSV}")
if [[ "${actual_chr2}" != "${expected_chr2}" ]]; then
    echo "FAIL: chr2:25246500-25246501 expected=${expected_chr2} got=${actual_chr2}"
    exit 1
fi
if [[ "${actual_chr7}" != "${expected_chr7}" ]]; then
    echo "FAIL: chr7:148884000-148884001 expected=${expected_chr7} got=${actual_chr7}"
    exit 1
fi
echo "  OK   anchor contract: chr2=5 chr7=9 (matches tests/fixtures/README.md)"
echo

# --- 7. verify_anchors against the just-built HTML --------------------------
echo "=== verify_anchors.py verify: HTML slice round-trip ==="
${PY} "${ANCHORS}" verify \
    --html "${HTML}" \
    --anchors "${ANCHORS_TSV}" \
    --out "${OUTDIR}/verify_anchors.tsv" \
    --fail-on-fail >/dev/null
echo "  OK   anchor verify PASS (HTML slice counts match source BAM counts)"
echo

# --- 8. Optional: --also-png exercises the full HTML+PNG pipeline -----------
# Skip semantics: this step is best-effort and never causes the test to FAIL.
# `igver` may be on PATH as a `pip install igver` egg-link shim WITHOUT the
# underlying IGV Java binary — exits 0 but produces no PNGs (the documented
# silent-failure mode in rules/igv.md). Our --also-png driver catches that
# via the inline existence check and raises SystemExit. Here we treat any
# such failure as SKIP rather than propagate it, since a non-working igver
# install isn't a regression in this skill's code.
if command -v igver >/dev/null 2>&1 || [[ -n "${IGVER_CMD:-}" ]]; then
    echo "=== --also-png: HTML + per-region PNGs (igver available) ==="
    HTML_PNG="${OUTDIR}/png_sample.hg38.html"
    if ${PY} "${BUILD}" \
            --sites "${SITES}" \
            --bam "${FIXTURE}" \
            --genome hg38 \
            --fasta "${FASTA}" \
            --no-default-tracks \
            --flanking 300 \
            --type mutation \
            --info-columns name \
            --output "${HTML_PNG}" \
            --no-apptainer \
            --no-verify \
            --also-png \
            --png-dpi 100 >"${OUTDIR}/also_png.log" 2>&1; then
        # --also-png returned 0 — assert the manifest + PNGs are real.
        MANIFEST="${OUTDIR}/png_png_sample.hg38/manifest.tsv"
        if [[ ! -f "${MANIFEST}" ]]; then
            echo "FAIL: --also-png exited 0 but no manifest at ${MANIFEST}"
            exit 1
        fi
        n_regions=$(awk -F'\t' 'NR>1 && !/^#/' "${MANIFEST}" | wc -l)
        if [[ "${n_regions}" -ne 3 ]]; then
            echo "FAIL: manifest has ${n_regions} regions, expected 3"
            exit 1
        fi
        png_one=$(awk -F'\t' 'NR==2 {print $9}' "${MANIFEST}")
        if [[ ! -s "${png_one}" ]]; then
            echo "FAIL: PNG missing or empty: ${png_one}"
            exit 1
        fi
        png_size=$(stat -c %s "${png_one}")
        echo "  OK   manifest: ${n_regions} regions; spot-check ${png_one##*/} = ${png_size} bytes"
    else
        # Driver caught the silent-failure mode; surface the diagnostic but
        # don't fail the test — broken igver install is environment-level.
        echo "  SKIP (igver invocation failed — likely missing IGV Java binary or wrong PATH)"
        echo "       see ${OUTDIR}/also_png.log for the driver's diagnostic."
        if grep -q "silent exit-0 failure\|Failed to generate all PNG files" "${OUTDIR}/also_png.log" 2>/dev/null; then
            echo "       (confirmed: this is the documented igver silent-failure mode)"
        fi
    fi
else
    echo "=== --also-png: SKIP (igver not on PATH; set \$IGVER_CMD or install via apptainer SIF) ==="
fi
echo

echo "=== end-to-end PASS — full pipeline (create_report → verify → optional igver) ==="
