#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/08_compartment_switch.sh
#
# Compartment A/B switches between two conditions. Run 06 first.
#
#   bash 08_compartment_switch.sh                       # day1 vs day56
#   COMPARE_LEFT=day7 COMPARE_RIGHT=day56 bash 08_compartment_switch.sh
#
# One library per condition is used, the first one listed for that condition in
# config/samples.tsv.
#
# Phasing
#   If 06 was run without PHASING_TRACK, point SWITCH_SIGNAL_DIR at a directory
#   holding <label>_signal.bed files (signal in the last column). Each sample is
#   then oriented so that compartment A has the higher signal, which is the rule
#   the earlier A/B assignment used.
#
# Outputs under results/hic/compartment/
#   switch_<left>_vs_<right>.pdf
#   switch_<left>_vs_<right>.bed
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

LEFT="${COMPARE_LEFT:-day1}"
RIGHT="${COMPARE_RIGHT:-day56}"
OUT="${OUT_DIR}/hic/compartment"
mkdir -p "${OUT}"

pick() {
    # first label of a condition
    samples hic | awk -F'\t' -v c="$1" '$1==c {print $3; exit}'
}

args=()
for cond in "${LEFT}" "${RIGHT}"; do
    label="$(pick "${cond}")"
    if [[ -z "${label}" ]]; then
        echo "ERROR: no Hi-C library for condition '${cond}' in config/samples.tsv" >&2
        exit 1
    fi
    eigs="${OUT}/${label}_eigs.tsv"
    if [[ ! -f "${eigs}" ]]; then
        echo "ERROR: ${eigs} not found. Run 06_compartment_eigs.sh first." >&2
        exit 1
    fi

    signal=""
    if [[ -n "${SWITCH_SIGNAL_DIR:-}" ]]; then
        cand="${SWITCH_SIGNAL_DIR}/${label}_signal.bed"
        if [[ -f "${cand}" ]]; then
            signal=",${cand}"
        else
            echo "[warn] ${cand} not found, using the eigenvector sign as is" >&2
        fi
    fi
    args+=(--sample "${cond}=${eigs}${signal}")
done

python3 "$(dirname "$0")/08_compartment_switch.py" \
    "${args[@]}" \
    --compare "${LEFT}=${RIGHT}" \
    --out "${OUT}/switch_${LEFT}_vs_${RIGHT}.pdf" \
    --out-bed "${OUT}/switch_${LEFT}_vs_${RIGHT}.bed"

echo "done -> ${OUT}/switch_${LEFT}_vs_${RIGHT}.pdf"
