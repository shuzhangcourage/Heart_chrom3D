#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/07_compartment_saddle.sh
#
# Saddle plot and compartment strength, for every library produced by
# 06_compartment_eigs.sh. Run 06 first.
#
#   RESOLUTION=100000 bash 07_compartment_saddle.sh
#
# Outputs under results/hic/compartment/
#   saddle_<label>.pdf / .png     one figure per library
#   saddle_all.pdf                all libraries side by side
#   compartment_strength.tsv      AA+BB over AB+BA per library
#   compartment_strength.pdf      bar chart of the above
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-100000}"
OUT="${OUT_DIR}/hic/compartment"
mkdir -p "${OUT}"

args=()
while IFS=$'\t' read -r cond rep label mark input; do
    mcool="${DATA_DIR}/${input}"
    eigs="${OUT}/${label}_eigs.tsv"
    exp="${OUT}/${label}_expected.tsv"

    if [[ ! -f "${mcool}" ]]; then
        echo "[skip] ${label}: ${mcool} not found" >&2
        continue
    fi
    if [[ ! -f "${eigs}" || ! -f "${exp}" ]]; then
        echo "[skip] ${label}: run 06_compartment_eigs.sh first" >&2
        continue
    fi

    if [[ "${mcool}" == *.mcool ]]; then
        cool="${mcool}::/resolutions/${RESOLUTION}"
    else
        cool="${mcool}"
    fi

    args+=(--sample "${label}=${cool},${eigs},${exp}")
done < <(samples hic)

if (( ${#args[@]} == 0 )); then
    echo "ERROR: nothing to plot." >&2
    echo "       either the matrices are missing, or 06 has not been run yet," >&2
    echo "       or RESOLUTION does not match the one used in 06." >&2
    exit 1
fi

python3 "$(dirname "$0")/07_compartment_saddle.py" \
    "${args[@]}" \
    --resolution "${RESOLUTION}" \
    --threads "${THREADS}" \
    --out "${OUT}/saddle" \
    --out-stacked "${OUT}/saddle_all.pdf" \
    --out-strength "${OUT}/compartment_strength.tsv"

echo "done -> ${OUT}/"
