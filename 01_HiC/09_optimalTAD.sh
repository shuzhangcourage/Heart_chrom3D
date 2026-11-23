#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/09_optimalTAD.sh
#
# TAD calling with optimalTAD.
#
#   RESOLUTION=10000 WMIN=40000 WMAX=40000 bash 09_optimalTAD.sh
#
# The --chipseq track is only used to pick the optimal window size, so either a
# H3K4me3 or an H3K27ac bigWig is acceptable ($OPTIMALTAD_CHIP in config.sh).
# In the published analysis min and max were both fixed at 40000, which skips
# the search and reproduces the value that was finally used.
#
# Input   results/hic/matrix/<label>_matrix_1kb.mcool
# Output  results/hic/tad/<label>.domains.tad
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-10000}"
WMIN="${WMIN:-40000}"
WMAX="${WMAX:-40000}"
OUT="${OUT_DIR}/hic/tad"
mkdir -p "${OUT}"

if [[ ! -f "${OPTIMALTAD_CHIP}" ]]; then
    echo "ERROR: --chipseq track not found: ${OPTIMALTAD_CHIP}" >&2
    echo "       set OPTIMALTAD_CHIP in config/config.sh" >&2
    exit 1
fi

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    mcool="${DATA_DIR}/${input}"
    if [[ ! -f "${mcool}" ]]; then
        echo "[skip] ${label}: ${mcool} not found" >&2
        continue
    fi
    found=$((found + 1))

    if [[ "${mcool}" == *.mcool ]]; then
        cool="${mcool}::/resolutions/${RESOLUTION}"
    else
        cool="${mcool}"
    fi

    echo "[${label}] optimalTAD, ${RESOLUTION} bp, window ${WMIN}-${WMAX}" >&2
    optimalTAD run \
        --hic "${cool}" \
        --output "${OUT}/work_${label}" \
        --chipseq "${OPTIMALTAD_CHIP}" \
        --np "${THREADS}" \
        --resolution "${RESOLUTION}" \
        --log2_hic \
        --log2_chip \
        --zscore_chip \
        --mammal \
        --window_size_min "${WMIN}" \
        --window_size_max "${WMAX}"

    tad="$(find "${OUT}/work_${label}" -name '*.domains.tad' | head -n 1 || true)"
    if [[ -z "${tad}" ]]; then
        echo "[warn] ${label}: no *.domains.tad produced, check the log above" >&2
        continue
    fi
    cp "${tad}" "${OUT}/${label}.domains.tad"
    echo "[${label}] -> ${OUT}/${label}.domains.tad" >&2
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no Hi-C matrices found under ${DATA_DIR}" >&2
    exit 1
fi

echo "done -> ${OUT}/"
