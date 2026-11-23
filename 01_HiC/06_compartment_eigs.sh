#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/06_compartment_eigs.sh
#
# Compartment A/B eigenvectors and the cis expected table, for every Hi-C
# library. Step 07 consumes both outputs.
#
#   RESOLUTION=50000 bash 06_compartment_eigs.sh
#   RESOLUTION=100000 PHASING_TRACK=/path/gc_100kb.bw bash 06_compartment_eigs.sh
#
# Phasing
#   eigs-cis returns eigenvectors with an arbitrary sign. The phasing track
#   (usually GC content) fixes the orientation so that A compartments have a
#   positive eigenvector. If PHASING_TRACK is unset the sign is left as is and
#   step 08 can re-orient it from a signal track.
#
# Outputs, per library, under results/hic/compartment/
#   <label>_eigs.tsv        eigenvectors
#   <label>_expected.tsv    cis expected, needed by the saddle plot
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-50000}"
VIEW="${VIEW_BED:-${REF_DIR}/mm10_viewdf.bed}"
OUT="${OUT_DIR}/hic/compartment"
mkdir -p "${OUT}"

if [[ ! -f "${VIEW}" ]]; then
    echo "ERROR: view file not found: ${VIEW}" >&2
    echo "       set VIEW_BED, or place mm10_viewdf.bed in \${REF_DIR}" >&2
    echo "       the file needs four columns: chrom, start, end, name" >&2
    exit 1
fi

phasing=()
if [[ -n "${PHASING_TRACK:-}" ]]; then
    phasing=(--phasing-track "${PHASING_TRACK}")
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

    echo "[${label}] balancing (cis only) at ${RESOLUTION} bp" >&2
    cooler balance --cis-only --force "${cool}"

    echo "[${label}] eigs-cis" >&2
    cooltools eigs-cis \
        --view "${VIEW}" \
        ${phasing[@]+"${phasing[@]}"} \
        -o "${OUT}/${label}_eigs.tsv" \
        "${cool}"

    echo "[${label}] expected-cis" >&2
    cooltools expected-cis \
        --view "${VIEW}" \
        -p "${THREADS}" \
        -o "${OUT}/${label}_expected.tsv" \
        "${cool}"
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no Hi-C matrices found under ${DATA_DIR}" >&2
    exit 1
fi

echo "done -> ${OUT}/"
