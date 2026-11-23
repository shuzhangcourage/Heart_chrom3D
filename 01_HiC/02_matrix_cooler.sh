#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/02_matrix_cooler.sh
#
# Turn pairs files into balanced multi-resolution coolers.
#
# Input   results/hic/pairs/<label>.pairs
# Output  results/hic/matrix/<label>_matrix_1kb.mcool
#         results/hic/matrix/<label>_matrix_1kb.cool    (1 kb, kept for cooler)
#
# Step order matters: the .cool is built at the base resolution, then zoomify
# adds the coarser resolutions and balances them in one pass.
#
# Resolutions
#   The zoom levels are listed explicitly because the analyses below ask for
#   5 kb, 10 kb, 50 kb and 100 kb. A geometric series (-r 1000N, i.e. 1, 2, 4,
#   8 ... kb) does not contain any of them, and cooltools would then fail with
#   KeyError: 5000. Override with ZOOM="5000,10000,..." if other resolutions are
#   needed.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RES="${RES:-1000}"
# every resolution used by the downstream steps must appear here
ZOOM="${ZOOM:-5000,10000,25000,50000,100000,250000,500000,1000000}"
OUT="${OUT_DIR}/hic/matrix"
PAIRS="${OUT_DIR}/hic/pairs"
mkdir -p "${OUT}"

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    pairs="${PAIRS}/${label}.pairs"
    if [[ ! -f "${pairs}" ]]; then
        echo "[skip] ${label}: ${pairs} not found" >&2
        continue
    fi
    found=$((found + 1))

    gz="${OUT}/${label}.pairs.gz"
    if [[ ! -f "${gz}" ]]; then
        echo "[${label}] bgzip + pairix" >&2
        bgzip -@ "${THREADS}" -c "${pairs}" > "${gz}"
        pairix -f "${gz}"
    fi

    echo "[${label}] cooler cload at ${RES} bp" >&2
    cooler cload pairix \
        -p "${THREADS}" \
        "${GENOME}:${RES}" \
        "${gz}" \
        "${OUT}/${label}_matrix_1kb.cool"

    echo "[${label}] zoomify + balance" >&2
    cooler zoomify \
        --balance \
        -p "${THREADS}" \
        -r "${RES},${ZOOM}" \
        -o "${OUT}/${label}_matrix_1kb.mcool" \
        "${OUT}/${label}_matrix_1kb.cool"
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no pairs files under ${PAIRS}. Run 01_mapping_pairtools.sh first." >&2
    exit 1
fi

echo "done -> ${OUT}/"
