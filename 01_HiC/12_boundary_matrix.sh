#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/12_boundary_matrix.sh
#
# Matrix around boundaries, relative values, and k-means clustering.
#
#   RESOLUTION=10000 BEFORE=200000 bash 12_boundary_matrix.sh
#
# The second derivative of the insulation score is used as the signal, so the
# heatmap shows where the boundary centre sits inside each class. Values are
# made relative to the centre bin before clustering, and the cluster labels are
# written to a BED file that steps 13 and 14 consume.
#
# Centre column
#   computeMatrix writes 6 metadata columns, then one column per bin. With B
#   bases before the reference point and a bin size of BS the centre column is
#   7 + B/BS. It is computed here, not hard-coded, so B and BS can be changed
#   without breaking the relative step.
#
# Input   results/hic/boundary/<label>_IS_d2.bw
#         results/hic/boundary/<label>_boundary.bed
# Output  results/hic/boundary/matrix/<label>_bmatrix.gz            absolute
#         results/hic/boundary/matrix/<label>_bmatrix.relative.gz   relative
#         results/hic/boundary/matrix/<label>_clusters.bed          k-means labels
#         results/hic/boundary/matrix/<label>_cluster_heatmap.pdf
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-10000}"
BEFORE="${BEFORE:-200000}"
AFTER="${AFTER:-$((BEFORE + RESOLUTION))}"
KMEANS="${KMEANS:-3}"
OUT="${OUT_DIR}/hic/boundary"
MAT="${OUT}/matrix"
mkdir -p "${MAT}"

if (( BEFORE % RESOLUTION != 0 )); then
    echo "ERROR: BEFORE (${BEFORE}) must be a multiple of RESOLUTION (${RESOLUTION})" >&2
    exit 1
fi
CENTRE=$(( 7 + BEFORE / RESOLUTION ))
echo "centre column = ${CENTRE} (6 metadata + ${BEFORE}/${RESOLUTION} bins)" >&2

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    bw="${OUT}/${label}_IS_d2.bw"
    bed="${OUT}/${label}_boundary.bed"
    if [[ ! -f "${bw}" || ! -f "${bed}" ]]; then
        echo "[skip] ${label}: need ${bw} and ${bed} (run 10 and 11 first)" >&2
        continue
    fi
    found=$((found + 1))

    mat="${MAT}/${label}_bmatrix.gz"

    echo "[${label}] computeMatrix reference-point" >&2
    computeMatrix reference-point \
        -S "${bw}" \
        -R "${bed}" \
        --referencePoint center \
        -b "${BEFORE}" -a "${AFTER}" \
        -bs "${RESOLUTION}" \
        --skipZeros \
        --missingDataAsZero \
        -p "${THREADS}" \
        -o "${mat}" \
        --outFileNameMatrix "${MAT}/${label}_bmatrix.tab" \
        --outFileSortedRegions "${MAT}/${label}_sorted_by_input.bed"

    echo "[${label}] subtract the centre bin" >&2
    awk -v mid="${CENTRE}" 'BEGIN{OFS="\t"}
    NR==1 { print; next }
    {
        n = NF
        printf "%s", $1
        for (i = 2; i <= 6; i++) printf "\t%s", $i
        for (i = 7; i <= n; i++) printf "\t%.6f", $i - $mid
        printf "\n"
    }' "${MAT}/${label}_bmatrix.tab" > "${MAT}/${label}_bmatrix.relative.tab"
    gzip -f "${MAT}/${label}_bmatrix.relative.tab"
    rel="${MAT}/${label}_bmatrix.relative.tab.gz"

    echo "[${label}] plotHeatmap with k-means (k=${KMEANS})" >&2
    plotHeatmap -m "${rel}" \
        -out "${MAT}/${label}_cluster_heatmap.pdf" \
        --colorMap seismic \
        --kmeans "${KMEANS}" \
        --heatmapHeight 14 \
        --yAxisLabel "Second derivative IS" \
        --xAxisLabel "Distance from boundary (kb)" \
        --zMin -0.6 --zMax 0.1 \
        --outFileSortedRegions "${MAT}/${label}_clusters.bed"

    echo "[${label}] -> ${MAT}/${label}_clusters.bed" >&2
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no boundary bigWig/BED found under ${OUT}" >&2
    exit 1
fi

echo "done -> ${MAT}/"
