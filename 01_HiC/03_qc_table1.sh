#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/03_qc_table1.sh
#
# Builds Table 1 from the pairtools statistics written by step 01.
#
#   bash 03_qc_table1.sh
#
# Output
#   results/hic/qc/table1.tsv
#   results/hic/qc/table1.md
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

OUT="${OUT_DIR}/hic/qc"
STATS="${OUT_DIR}/hic/pairs"
mkdir -p "${OUT}"

args=()
while IFS=$'\t' read -r cond rep label mark input; do
    f="${STATS}/${label}.dedup.stats.txt"
    if [[ ! -f "${f}" ]]; then
        f="${STATS}/${label}.parse.stats.txt"
    fi
    if [[ ! -f "${f}" ]]; then
        echo "[skip] ${label}: no pairtools stats found" >&2
        continue
    fi
    args+=(--sample "${label}=${f}")
done < <(samples hic)

if (( ${#args[@]} == 0 )); then
    echo "ERROR: no statistics files under ${STATS}" >&2
    echo "       run 01_mapping_pairtools.sh first" >&2
    exit 1
fi

python3 "$(dirname "$0")/03_qc_table1.py" \
    "${args[@]}" \
    --out "${OUT}/table1.tsv" \
    --out-markdown "${OUT}/table1.md"

echo "done -> ${OUT}/table1.tsv"
