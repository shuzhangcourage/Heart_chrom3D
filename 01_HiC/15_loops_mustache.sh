#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/15_loops_mustache.sh
#
# Chromatin loop calling with mustache, autosomes only.
#
#   RESOLUTION=5kb PT=0.01 bash 15_loops_mustache.sh
#
# Input   results/hic/matrix/<label>_matrix_1kb.mcool
# Output  results/hic/loops/<label>_loops.tsv
#
# Parameters from the manuscript: -r 5kb, -pt 0.01, chromosomes chr1-chr19.
# Restricting to autosomes is deliberate and excludes chrX and chrY.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-5kb}"
PT="${PT:-${MUSTACHE_PT:-0.01}}"
OUT="${OUT_DIR}/hic/loops"
mkdir -p "${OUT}"

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19)

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    mcool="${DATA_DIR}/${input}"
    if [[ ! -f "${mcool}" ]]; then
        echo "[skip] ${label}: ${mcool} not found" >&2
        continue
    fi
    found=$((found + 1))

    echo "[${label}] mustache, ${RESOLUTION}, pt ${PT}" >&2
    mustache \
        -f "${mcool}" \
        -ch "${CHROMS[@]}" \
        -r "${RESOLUTION}" \
        -pt "${PT}" \
        -p "${THREADS}" \
        -o "${OUT}/${label}_loops.tsv"
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no Hi-C matrices found under ${DATA_DIR}" >&2
    exit 1
fi

echo "done -> ${OUT}/"
