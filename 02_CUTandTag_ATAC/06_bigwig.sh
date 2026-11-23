#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02_CUTandTag_ATAC/06_bigwig.sh
#
# Coverage tracks for visualisation and for the correlation analysis in step 09.
#
#   bash 06_bigwig.sh
#
# Input   results/cuttag/bam/<label>.filt.bam
# Output  results/cuttag/bigwig/<label>.bw
#
# CPM normalisation is used, so tracks are comparable within a library but not
# necessarily across libraries with very different signal-to-background. If
# spike-in factors are available (step 08), re-run with
#
#   SPIKEIN_TABLE=results/cuttag/spikein_factors.tsv bash 06_bigwig.sh
#
# and each track is additionally multiplied by its spike-in factor, writing to
# results/cuttag/bigwig_spikein/<label>.bw instead.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

BAM="${OUT_DIR}/cuttag/bam"
SPIKEIN_TABLE="${SPIKEIN_TABLE:-}"

if [[ -n "${SPIKEIN_TABLE}" && -f "${SPIKEIN_TABLE}" ]]; then
    BW="${OUT_DIR}/cuttag/bigwig_spikein"
    echo "spike-in scaling: ${SPIKEIN_TABLE}" >&2
else
    BW="${OUT_DIR}/cuttag/bigwig"
    SPIKEIN_TABLE=""
fi
mkdir -p "${BW}"

lookup_factor() {
    # column 2 of the factors table for one label, or empty
    awk -F'\t' -v l="$1" '$1==l {print $2; exit}' "${SPIKEIN_TABLE}"
}

found=0
for dt in cuttag atac; do
    while IFS=$'\t' read -r cond rep label mark input; do
        in_bam="${BAM}/${label}.filt.bam"
        if [[ ! -f "${in_bam}" ]]; then
            echo "[skip] ${label}: ${in_bam} not found (run 05 first)" >&2
            continue
        fi
        found=$((found + 1))

        extra=()
        if [[ -n "${SPIKEIN_TABLE}" ]]; then
            factor="$(lookup_factor "${label}")"
            if [[ -z "${factor}" ]]; then
                echo "[skip] ${label}: no spike-in factor, not in ${SPIKEIN_TABLE}" >&2
                continue
            fi
            extra=(--scaleFactor "${factor}")
        fi

        echo "[${label}] bamCoverage -> ${BW}/${label}.bw" >&2
        bamCoverage \
            -b "${in_bam}" \
            -o "${BW}/${label}.bw" \
            --binSize 10 \
            --normalizeUsing CPM \
            --extendReads \
            --numberOfProcessors "${THREADS}" \
            ${extra[@]+"${extra[@]}"}
    done < <(samples "${dt}")
done

if (( found == 0 )); then
    echo "ERROR: no filtered BAM under ${BAM}. Run 05_filter.sh first." >&2
    exit 1
fi

echo "done -> ${BW}/"
