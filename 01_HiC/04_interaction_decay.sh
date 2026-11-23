#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/04_interaction_decay.sh
#
# Intra-chromosomal interaction decay for every Hi-C library in
# config/samples.tsv.
#
# Legend labels are <cond> for merged libraries and <cond>_rep<N> otherwise, so
# replicates stay distinguishable. To reproduce a figure with one curve per
# condition, delete the unused rows from config/samples.tsv.
#
#   RESOLUTION=5000 VIEW_BED=/path/to/mm10_viewdf.bed bash 04_interaction_decay.sh
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-5000}"
OUT="${OUT_DIR}/hic"
mkdir -p "${OUT}"

args=()
while IFS=$'\t' read -r cond rep label mark input; do
    f="${DATA_DIR}/${input}"
    if [[ ! -f "${f}" ]]; then
        echo "[skip] ${label}: ${f} not found" >&2
        continue
    fi
    if [[ "${rep}" == "merged" ]]; then
        name="${cond}"
    else
        name="${cond}_rep${rep}"
    fi
    args+=(--sample "${name}=${f}")
done < <(samples hic)

if (( ${#args[@]} == 0 )); then
    echo "ERROR: no Hi-C matrices found under ${DATA_DIR}" >&2
    echo "       check config/config.sh and config/samples.tsv" >&2
    exit 1
fi

extra=()
if [[ -n "${VIEW_BED:-}" ]]; then
    extra+=(--view "${VIEW_BED}")
fi

python3 "$(dirname "$0")/04_interaction_decay.py" \
    "${args[@]}" \
    --resolution "${RESOLUTION}" \
    --threads "${THREADS}" \
    "${extra[@]}" \
    --out "${OUT}/interaction_decay.pdf" \
    --out-table "${OUT}/interaction_decay.tsv"

echo "done -> ${OUT}/interaction_decay.pdf"
