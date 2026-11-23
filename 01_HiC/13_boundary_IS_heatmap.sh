#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/13_boundary_IS_heatmap.sh
#
# Insulation score heatmap per boundary class, keeping the order that step 12
# produced. Two figures per library: absolute insulation score, and insulation
# score made relative to the centre bin.
#
#   RESOLUTION=10000 BEFORE=200000 bash 13_boundary_IS_heatmap.sh
#
# Input   results/hic/boundary/<label>_IS.bw
#         results/hic/boundary/matrix/<label>_clusters.bed
# Output  results/hic/boundary/matrix/<label>_IS_sorted_heatmap.pdf
#         results/hic/boundary/matrix/<label>_IS_sorted_relative_heatmap.pdf
#
# The centre column is computed from BEFORE and RESOLUTION, as in step 12.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-10000}"
BEFORE="${BEFORE:-200000}"
AFTER="${AFTER:-$((BEFORE + RESOLUTION))}"
OUT="${OUT_DIR}/hic/boundary"
MAT="${OUT}/matrix"
mkdir -p "${MAT}"

if (( BEFORE % RESOLUTION != 0 )); then
    echo "ERROR: BEFORE (${BEFORE}) must be a multiple of RESOLUTION (${RESOLUTION})" >&2
    exit 1
fi
CENTRE=$(( 7 + BEFORE / RESOLUTION ))

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    bw="${OUT}/${label}_IS.bw"
    regions="${MAT}/${label}_clusters.bed"
    if [[ ! -f "${bw}" || ! -f "${regions}" ]]; then
        echo "[skip] ${label}: need ${bw} and ${regions} (run 10 and 12 first)" >&2
        continue
    fi
    found=$((found + 1))

    mat="${MAT}/${label}_IS_sorted.gz"

    echo "[${label}] computeMatrix on the clustered regions" >&2
    computeMatrix reference-point \
        -S "${bw}" \
        -R "${regions}" \
        --referencePoint center \
        -b "${BEFORE}" -a "${AFTER}" \
        -bs "${RESOLUTION}" \
        --skipZeros \
        --missingDataAsZero \
        -p "${THREADS}" \
        -o "${mat}" \
        --outFileNameMatrix "${MAT}/${label}_IS_sorted.tab" \
        --outFileSortedRegions "${MAT}/${label}_IS_sorted_regions.bed"

    echo "[${label}] heatmap, absolute" >&2
    plotHeatmap -m "${mat}" \
        -out "${MAT}/${label}_IS_sorted_heatmap.pdf" \
        --colorMap RdBu_r \
        --heatmapHeight 14 \
        --sortRegions no \
        --whatToShow 'plot, heatmap and colorbar' \
        --yAxisLabel "Insulation score" \
        --samplesLabel IS \
        --yMin -1.0 --yMax 0.5 \
        --zMin -0.6 --zMax 0.2 \
        --legendLocation none

    echo "[${label}] heatmap, relative to the centre bin" >&2
    awk -v mid="${CENTRE}" 'BEGIN{OFS="\t"}
    NR==1 { print; next }
    {
        n = NF
        printf "%s", $1
        for (i = 2; i <= 6; i++) printf "\t%s", $i
        for (i = 7; i <= n; i++) printf "\t%.6f", $i - $mid
        printf "\n"
    }' "${MAT}/${label}_IS_sorted.tab" > "${MAT}/${label}_IS_sorted_relative.tab"
    gzip -f "${MAT}/${label}_IS_sorted_relative.tab"

    plotHeatmap -m "${MAT}/${label}_IS_sorted_relative.tab.gz" \
        -out "${MAT}/${label}_IS_sorted_relative_heatmap.pdf" \
        --colorMap RdBu_r \
        --heatmapHeight 14 \
        --sortRegions no \
        --whatToShow 'plot, heatmap and colorbar' \
        --yAxisLabel "relative IS" \
        --xAxisLabel "Distance from boundary (kb)" \
        --samplesLabel "relative IS" \
        --yMin -0.2 --yMax 1.0 \
        --zMin -0.4 --zMax 0.8 \
        --legendLocation none
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: nothing to plot. Run 10_insulation_IS.sh and 12_boundary_matrix.sh first." >&2
    exit 1
fi

echo "done -> ${MAT}/"
