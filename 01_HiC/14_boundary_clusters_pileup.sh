#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/14_boundary_clusters_pileup.sh
#
# Pile-up of observed/expected contacts around the boundaries of each class.
#
#   RESOLUTION=10000 FLANK=500000 bash 14_boundary_clusters_pileup.sh
#
# The cluster labels from step 12 are split into one BED per class, then
# coolpuppy accumulates local contacts for each class.
#
# Input   results/hic/boundary/matrix/<label>_clusters.bed
#         results/hic/matrix/<label>_matrix_1kb.mcool
# Output  results/hic/pileup/<label>_cluster_<N>.bed
#         results/hic/pileup/<label>_cluster_<N>_pileup.pdf
#         results/hic/pileup/<label>_expected_<RES>kb.tsv   (if not present)
#
# Requires coolpup.py and plotpup.py on PATH.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

RESOLUTION="${RESOLUTION:-10000}"
FLANK="${FLANK:-500000}"
OUT="${OUT_DIR}/hic/pileup"
MAT="${OUT_DIR}/hic/boundary/matrix"
mkdir -p "${OUT}"

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    clusters="${MAT}/${label}_clusters.bed"
    mcool="${DATA_DIR}/${input}"
    if [[ ! -f "${clusters}" ]]; then
        echo "[skip] ${label}: ${clusters} not found (run 12 first)" >&2
        continue
    fi
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

    expected="${OUT}/${label}_expected_$((RESOLUTION / 1000))kb.tsv"
    if [[ ! -f "${expected}" ]]; then
        echo "[${label}] expected-cis at ${RESOLUTION} bp" >&2
        cooltools expected-cis -p "${THREADS}" -o "${expected}" "${cool}"
    fi

    # one BED per cluster; the cluster label is the last column
    # remove last run's files first, the awk below appends
    rm -f "${OUT}/${label}"_cluster_*.bed
    awk -v out="${OUT}/${label}" 'BEGIN{OFS="\t"}
    $NF ~ /^cluster_/ {
        f = out "_" $NF ".bed"
        print $1, $2, $3 >> f
    }' "${clusters}"

    for bed in "${OUT}/${label}"_cluster*.bed; do
        [[ -f "${bed}" ]] || continue
        prefix="$(basename "${bed}" .bed)"

        echo "[${label}] pile-up ${prefix}" >&2
        coolpup.py "${cool}" "${bed}" \
            --features_format bed \
            --local \
            --nshifts 0 \
            --ignore_diags 0 \
            --flank "${FLANK}" \
            --expected "${expected}" \
            --outname "${OUT}/${prefix}_pileup.clpy" \
            --nproc "${THREADS}"

        plotpup.py \
            --cmap gist_ncar \
            --no_score \
            --not_symmetrical \
            --dpi 600 \
            --height 2 \
            --vmax 1.4 --vmin 0.6 \
            --input_pups "${OUT}/${prefix}_pileup.clpy" \
            --output "${OUT}/${prefix}_pileup.pdf"
    done
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: nothing to pile up. Run 12_boundary_matrix.sh first." >&2
    exit 1
fi

echo "done -> ${OUT}/"
