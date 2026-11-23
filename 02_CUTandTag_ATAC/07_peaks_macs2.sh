#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02_CUTandTag_ATAC/07_peaks_macs2.sh
#
# Peak calling with MACS2. The significance cutoff is a q-value, set by
# MACS2_QVALUE in config/config.sh (default 0.05).
#
#   bash 07_peaks_macs2.sh
#
# Input   results/cuttag/bam/<label>.filt.bam
#         optional control BAM, from the "input" column of config/samples.tsv
# Output  results/cuttag/peaks/<label>_peaks.narrowPeak
#         results/cuttag/peaks/<label>_summits.bed
#
# -f BAMPE lets MACS2 use the actual fragment length from the paired reads
# instead of building a model from the read ends. An IgG or input library is
# used as control only when the samples table names one for that library.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

BAM="${OUT_DIR}/cuttag/bam"
PEAK="${OUT_DIR}/cuttag/peaks"
QVALUE="${MACS2_QVALUE:-0.05}"
MACS2_EXTRA="${MACS2_EXTRA:-}"
mkdir -p "${PEAK}"

found=0
for dt in cuttag atac; do
    while IFS=$'\t' read -r cond rep label mark input; do
        in_bam="${BAM}/${label}.filt.bam"
        if [[ ! -f "${in_bam}" ]]; then
            echo "[skip] ${label}: ${in_bam} not found (run 05 first)" >&2
            continue
        fi
        found=$((found + 1))

        control=()
        if [[ -n "${input}" && "${input}" != "-" ]]; then
            ctrl_bam="${DATA_DIR}/${input}"
            if [[ -f "${ctrl_bam}" ]]; then
                control=(-c "${ctrl_bam}")
            else
                echo "[${label}] control ${ctrl_bam} not found, calling without it" >&2
            fi
        fi

        echo "[${label}] macs2 callpeak -q ${QVALUE}" >&2
        macs2 callpeak \
            -t "${in_bam}" \
            ${control[@]+"${control[@]}"} \
            -f BAMPE \
            -g mm \
            -q "${QVALUE}" \
            --keep-dup all \
            -n "${label}" \
            --outdir "${PEAK}" \
            ${MACS2_EXTRA}
    done < <(samples "${dt}")
done

if (( found == 0 )); then
    echo "ERROR: no filtered BAM under ${BAM}. Run 05_filter.sh first." >&2
    exit 1
fi

echo "done -> ${PEAK}/"
