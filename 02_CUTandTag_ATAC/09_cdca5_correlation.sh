#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02_CUTandTag_ATAC/09_cdca5_correlation.sh
#
# Correlation between Cdca5 occupancy and the chromatin marks, in 100 kb bins.
#
#   bash 09_cdca5_correlation.sh
#
# Input   results/cuttag/bigwig/<label>.bw            from step 06
#         public ChIP-seq bigWigs listed in config/samples.tsv (source=public)
# Output  results/correlation/counts_100kb.npz        binned signal
#         results/correlation/correlation.tsv         Spearman matrix
#         results/correlation/correlation_matrix.pdf  all tracks, all-vs-all
#         results/correlation/cdca5_correlation.pdf   Cdca5 vs marks, per condition
#
# This is the only step that mixes data types: CUT&Tag, ATAC-seq and the public
# histone ChIP-seq tracks are binned on the same grid so that a single
# correlation matrix can be computed.
#
# BIN_SIZE is 100 kb. The public tracks are already coverage bigWigs, so no
# further normalisation is applied here; multiBigwigSummary averages signal
# inside each bin.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

BIN_SIZE="${BIN_SIZE:-100000}"
COR="${OUT_DIR}/correlation"
mkdir -p "${COR}"

BWS=()
LABELS=()
META="${COR}/tracks.tsv"
N_TRACKS=0
printf 'label\tcond\tmark\n' > "${META}"

add_track() {
    # add one bigWig to the list if it exists
    local bw="$1" label="$2" cond="$3" mark="$4"
    if [[ ! -f "${bw}" ]]; then
        echo "[skip] ${label}: ${bw} not found" >&2
        return
    fi
    BWS+=("${bw}")
    LABELS+=("${label}")
    N_TRACKS=$((N_TRACKS + 1))
    printf '%s\t%s\t%s\n' "${label}" "${cond}" "${mark}" >> "${META}"
}

# --- our own CUT&Tag and ATAC-seq tracks ------------------------------------
for dt in cuttag atac; do
    while IFS=$'\t' read -r cond rep label mark input; do
        add_track "${OUT_DIR}/cuttag/bigwig/${label}.bw" "${label}" "${cond}" "${mark}"
    done < <(samples "${dt}")
done

# --- public histone ChIP-seq tracks ----------------------------------------
while IFS=$'\t' read -r cond rep label mark input; do
    add_track "${DATA_DIR}/${input}" "${label}" "${cond}" "${mark}"
done < <(samples chipseq)

if (( N_TRACKS < 2 )); then
    echo "ERROR: fewer than two bigWigs found; run 06_bigwig.sh first." >&2
    exit 1
fi

echo "summarising ${N_TRACKS} tracks in ${BIN_SIZE} bp bins" >&2
multiBigwigSummary bins \
    -b "${BWS[@]}" \
    --labels "${LABELS[@]}" \
    --binSize "${BIN_SIZE}" \
    --numberOfProcessors "${THREADS}" \
    -o "${COR}/counts_100kb.npz"

python3 "$(dirname "$0")/09_cdca5_correlation.py" \
    --matrix "${COR}/counts_100kb.npz" \
    --meta "${META}" \
    --outdir "${COR}"

echo "done -> ${COR}/"
