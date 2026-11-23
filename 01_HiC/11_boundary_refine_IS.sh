#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/11_boundary_refine_IS.sh
#
# Refine TAD edges into boundary bins. Run 09 and 10 first.
#
#   RESOLUTION=10000 bash 11_boundary_refine_IS.sh
#
# Input   results/hic/boundary/<label>_IS.tsv
#         results/hic/tad/<label>.domains.tad
# Output  results/hic/boundary/<label>_boundary.tsv
#         results/hic/boundary/<label>_boundary.bed
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-10000}"
WINDOW="${WINDOW:-40000}"
OUT="${OUT_DIR}/hic/boundary"
TAD="${OUT_DIR}/hic/tad"
mkdir -p "${OUT}"

extra=()
if [[ "${DROP_POSITIVE:-0}" == "1" ]]; then
    extra+=(--drop-positive)
fi

args=()
while IFS=$'\t' read -r cond rep label mark input; do
    is_tsv="${OUT}/${label}_IS.tsv"
    tad="${TAD}/${label}.domains.tad"
    if [[ ! -f "${is_tsv}" || ! -f "${tad}" ]]; then
        echo "[skip] ${label}: need ${is_tsv} and ${tad}" >&2
        continue
    fi
    args+=(--sample "${label}=${is_tsv}:${tad}")
done < <(samples hic)

if (( ${#args[@]} == 0 )); then
    echo "ERROR: nothing to refine. Run 09_optimalTAD.sh and 10_insulation_IS.sh first." >&2
    exit 1
fi

python3 "$(dirname "$0")/11_boundary_refine_IS.py" \
    "${args[@]}" \
    --bin-size "${RESOLUTION}" \
    --window "${WINDOW}" \
    ${extra[@]+"${extra[@]}"} \
    --out-dir "${OUT}"

echo "done -> ${OUT}/"
