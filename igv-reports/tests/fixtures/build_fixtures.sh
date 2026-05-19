#!/usr/bin/env bash
# build_fixtures.sh — regenerate tests/fixtures/tiny_colo829.hg38.bam from
# the publicly released ONT COLO829BL reads.
#
# The output BAM is committed to the repo (it's small public data — see
# fixtures/README.md). Regenerate only when you need to expand the slice
# regions, change subsample rate, or update for a new basecaller version.
# If the output counts change, also update tests/smoke/test_slice_count.py
# anchor constants and any integration scenarios.sh expected values.
set -euo pipefail

FIX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Resolve samtools.
if [[ -n "${SAMTOOLS_SIF:-}" && -f "${SAMTOOLS_SIF}" ]]; then
    BINDS=()
    if [[ -n "${IGV_REPORTS_BIND:-}" ]]; then
        IFS=':' read -ra BP <<< "${IGV_REPORTS_BIND}"
        for p in "${BP[@]}"; do [[ -d "$p" ]] && BINDS+=(--bind "$p"); done
    fi
    SAM=(apptainer exec --cleanenv "${BINDS[@]}" "${SAMTOOLS_SIF}" samtools)
elif command -v samtools >/dev/null 2>&1; then
    SAM=(samtools)
else
    echo "ERROR: no samtools available — install (pip/conda) or set SAMTOOLS_SIF" >&2
    exit 1
fi

# Source BAM — must be supplied via env. Public source: ENA project PRJEB57425
# (ONT COLO829 release).
SRC="${COLO829BL_BAM:-}"
if [[ -z "${SRC}" || ! -f "${SRC}" ]]; then
    echo "ERROR: source BAM not provided." >&2
    echo "       Set COLO829BL_BAM=<path-to-COLO829BL-ONT-BAM> and re-run." >&2
    echo "       Public source: ENA project PRJEB57425." >&2
    exit 1
fi

OUT="${FIX_DIR}/tiny_colo829.hg38.bam"

echo "[build_fixtures] source: ${SRC}"
echo "[build_fixtures] output: ${OUT}"
echo "[build_fixtures] regions: chr2:25245000-25248000 (DNMT3A), chr7:148882000-148886000 (EZH2)"
echo "[build_fixtures] subsample: 0.2, seed 42"

"${SAM[@]}" view -bh -F 1536 --subsample 0.2 --subsample-seed 42 \
    "${SRC}" \
    chr2:25245000-25248000 chr7:148882000-148886000 \
    -o "${OUT}"
"${SAM[@]}" index "${OUT}"

echo "[build_fixtures] sizes:"
ls -lh "${OUT}" "${OUT}.bai"

echo "[build_fixtures] anchor counts (must remain stable across regens):"
chr2_n=$("${SAM[@]}" view -c -F 1536 "${OUT}" chr2:25246500-25246501)
chr7_n=$("${SAM[@]}" view -c -F 1536 "${OUT}" chr7:148884000-148884001)
echo "  chr2:25246500-25246501 = ${chr2_n}"
echo "  chr7:148884000-148884001 = ${chr7_n}"

if [[ "${chr2_n}" != "5" || "${chr7_n}" != "9" ]]; then
    echo
    echo "WARNING: anchor counts have changed from the committed fixture's contract" >&2
    echo "         (chr2=5, chr7=9). Update tests/smoke/test_slice_count.py and any" >&2
    echo "         integration scenarios.sh expected values, then commit both the new" >&2
    echo "         BAM and the updated test constants together." >&2
fi
