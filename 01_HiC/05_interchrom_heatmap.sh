#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/05_interchrom_heatmap.sh
#
# Inter-chromosome contact heatmap, one figure per Hi-C library.
#
#   RESOLUTION=5000 bash 05_interchrom_heatmap.sh
#   CHROMS=chr1,chr2 bash 05_interchrom_heatmap.sh     # subset for a quick look
#
# The published figure is drawn at 5 kb. Pass RESOLUTION=500000 for a 500 kb
# rendering of the same matrix.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-5000}"
OUT="${OUT_DIR}/hic/interchrom"
mkdir -p "${OUT}"

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    f="${DATA_DIR}/${input}"
    if [[ ! -f "${f}" ]]; then
        echo "[skip] ${label}: ${f} not found" >&2
        continue
    fi
    found=$((found + 1))

    chroms_arg=()
    if [[ -n "${CHROMS:-}" ]]; then
        chroms_arg=(--chroms "${CHROMS}")
    fi

    python3 "$(dirname "$0")/05_interchrom_heatmap.py" \
        --sample "${label}=${f}" \
        --resolution "${RESOLUTION}" \
        --threads "${THREADS}" \
        "${chroms_arg[@]}" \
        --out "${OUT}/interchrom_${label}.pdf" \
        --out-table "${OUT}/interchrom_${label}.tsv"
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no Hi-C matrices found under ${DATA_DIR}" >&2
    echo "       check config/config.sh and config/samples.tsv" >&2
    exit 1
fi

echo "done -> ${OUT}/"
