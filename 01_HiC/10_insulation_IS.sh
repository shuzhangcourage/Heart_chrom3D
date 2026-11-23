#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/10_insulation_IS.sh
#
# Insulation score and its first and second derivatives.
#
#   RESOLUTION=10000 WINDOW=40000 bash 10_insulation_IS.sh
#
# The insulation score is computed with cooltools at a fixed window, then
# converted to a bedGraph so the derivatives can be taken bin by bin. The
# second derivative peaks at the boundary centre, which is what step 12 uses to
# place the matrix window.
#
# Input   results/hic/matrix/<label>_matrix_1kb.mcool
# Output  results/hic/boundary/<label>_IS.tsv          cooltools table
#         results/hic/boundary/<label>_IS.bw           insulation score
#         results/hic/boundary/<label>_IS_d1.bw        first derivative
#         results/hic/boundary/<label>_IS_d2.bw        second derivative
#
# Requires the UCSC tools bigWigToBedGraph and bedGraphToBigWig on PATH.
# Intermediate bedGraph files are removed unless KEEP_TMP=1.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-10000}"
WINDOW="${WINDOW:-40000}"
OUT="${OUT_DIR}/hic/boundary"
mkdir -p "${OUT}"

derivatives() {
    local bw="$1" prefix="$2"
    local bg="${prefix}.bg"
    local bg1="${prefix}_d1.bg"
    local bg2="${prefix}_d2.bg"

    bigWigToBedGraph "${bw}" "${bg}"

    # first derivative: value[i] - value[i-1] within a chromosome
    awk 'BEGIN{OFS="\t"}
    {
        if (NR==1) { pc=$1; ps=$2; pe=$3; pv=$4; next }
        if ($1==pc) print pc, ps, pe, $4-pv
        pc=$1; ps=$2; pe=$3; pv=$4
    }' "${bg}" | sort -k1,1 -k2,2n > "${bg1}"

    # second derivative: value[i] - 2*value[i-1] + value[i-2]
    awk 'BEGIN{OFS="\t"}
    {
        if (NR==1) { c1=$1; s1=$2; e1=$3; v1=$4; next }
        if (NR==2) { c2=$1; s2=$2; e2=$3; v2=$4; next }
        if ($1==c1 && c1==c2) print c2, s2, e2, $4-2*v2+v1
        c1=c2; s1=s2; e1=e2; v1=v2
        c2=$1; s2=$2; e2=$3; v2=$4
    }' "${bg}" | sort -k1,1 -k2,2n > "${bg2}"

    bedGraphToBigWig "${bg1}" "${GENOME}" "${prefix}_d1.bw"
    bedGraphToBigWig "${bg2}" "${GENOME}" "${prefix}_d2.bw"

    if [[ "${KEEP_TMP:-0}" != "1" ]]; then
        rm -f "${bg}" "${bg1}" "${bg2}"
    fi
}

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

    echo "[${label}] cooltools insulation, ${WINDOW} bp window" >&2
    cooltools insulation \
        "${cool}" \
        "${WINDOW}" \
        -p "${THREADS}" \
        -o "${OUT}/${label}_IS.tsv" \
        --bigwig

    bw="${OUT}/${label}_IS.bw"
    if [[ ! -f "${bw}" ]]; then
        # cooltools names the bigWig after the -o prefix
        alt="${OUT}/${label}_IS_${WINDOW}.bw"
        if [[ -f "${alt}" ]]; then
            mv "${alt}" "${bw}"
        else
            echo "[warn] ${label}: bigWig not found, skipping derivatives" >&2
            continue
        fi
    fi

    echo "[${label}] derivatives" >&2
    derivatives "${bw}" "${OUT}/${label}_IS"
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no Hi-C matrices found under ${DATA_DIR}" >&2
    exit 1
fi

echo "done -> ${OUT}/"
