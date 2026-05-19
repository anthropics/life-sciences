#!/usr/bin/env bash
# tests/run_all.sh — orchestrate the three test layers in order.
#
# Author: Samuel Ahuno
# Purpose:
#   1. unit (~1 s)      — pure-Python parser tests; pytest.
#   2. smoke (~3 s)     — samtools subprocess + slice-decode round-trip
#                          against the committed fixture; pytest.
#   3. integration      — full cohort build + verify-cohort + verify-anchors
#                          end-to-end; bash scenarios.sh under each demo.
#                          Skipped (exit 77) when the IGV_REPORTS_TEST_BAM_*
#                          env vars are unset AND the MSKCC default paths
#                          don't exist.
#
# Usage:
#   bash tests/run_all.sh              # all three layers
#   bash tests/run_all.sh --unit-only  # layer 1 only — instant feedback
#   bash tests/run_all.sh --no-integration  # layers 1 + 2 (fast everywhere)
#   bash tests/run_all.sh --integration-only  # layer 3 only — for the slow lane
#
# Exit code:
#   0  — every requested layer passed (or was legitimately skipped).
#   1+ — at least one layer failed; output preserved for debugging.
set -euo pipefail

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_DIR="$(cd "${TESTS_DIR}/.." && pwd)"

RUN_UNIT=1
RUN_SMOKE=1
RUN_INTEGRATION=1

for arg in "$@"; do
    case "$arg" in
        --unit-only)         RUN_SMOKE=0; RUN_INTEGRATION=0 ;;
        --no-integration)    RUN_INTEGRATION=0 ;;
        --integration-only)  RUN_UNIT=0; RUN_SMOKE=0 ;;
        -h|--help)
            sed -n '3,20p' "$0"
            exit 0
            ;;
        *)
            echo "ERROR: unknown flag: $arg" >&2
            echo "       Use --help to see options." >&2
            exit 2
            ;;
    esac
done

# Pick a Python with pytest. Prefer the snakemake conda env (where all
# project tooling lives); fall back to PATH `python3`.
PY="${IGV_REPORTS_PY:-}"
if [[ -z "${PY}" ]]; then
    if [[ -x /home/ahunos/miniforge3/envs/snakemake/bin/python ]]; then
        PY=/home/ahunos/miniforge3/envs/snakemake/bin/python
    elif command -v python3 >/dev/null 2>&1; then
        PY=$(command -v python3)
    else
        echo "ERROR: no python3 available. Set IGV_REPORTS_PY=<path-to-python>" >&2
        exit 2
    fi
fi

FAILS=0
SKIPS=0

run_layer() {
    local name="$1"; shift
    local desc="$1"; shift
    echo "=== ${name}: ${desc} ==="
    if "$@"; then
        echo "    ${name} PASS"
    else
        local rc=$?
        if [[ $rc -eq 77 ]]; then
            echo "    ${name} SKIP (exit 77 — see message above)"
            SKIPS=$((SKIPS + 1))
        else
            echo "    ${name} FAIL (exit ${rc})"
            FAILS=$((FAILS + 1))
        fi
    fi
    echo
}

# --- Layer 1: unit ---------------------------------------------------------
if [[ $RUN_UNIT -eq 1 ]]; then
    run_layer "unit" "pure-Python parsers" \
        "${PY}" -m pytest "${TESTS_DIR}/unit/" -q
fi

# --- Layer 2: smoke --------------------------------------------------------
if [[ $RUN_SMOKE -eq 1 ]]; then
    run_layer "smoke" "samtools + slice-decode round-trip" \
        "${PY}" -m pytest "${TESTS_DIR}/smoke/" -q
fi

# --- Layer 3: integration --------------------------------------------------
# Each scenarios.sh exits 77 if its required BAMs aren't available; we treat
# that as a skip rather than a failure so the suite is portable.
if [[ $RUN_INTEGRATION -eq 1 ]]; then
    # end_to_end: uses the committed tiny_colo829 fixture (~30 s, runs in CI).
    run_layer "integration / end_to_end" "full pipeline against committed fixture" \
        bash "${TESTS_DIR}/integration/end_to_end/scenarios.sh"
    run_layer "integration / cohort_verify" "cohort structural verifier scenarios" \
        bash "${TESTS_DIR}/integration/cohort_verify/scenarios.sh"
    run_layer "integration / anchor_verify" "anchor content verifier scenarios" \
        bash "${TESTS_DIR}/integration/anchor_verify/scenarios.sh"
fi

echo "=== summary ==="
echo "  failures: ${FAILS}"
echo "  skips:    ${SKIPS}"

if [[ $FAILS -gt 0 ]]; then
    exit 1
fi
exit 0
