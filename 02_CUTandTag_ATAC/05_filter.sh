#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02_CUTandTag_ATAC/05_filter.sh
#
# Mapping-quality filter. Everything downstream (coverage tracks, peak calling)
# reads this BAM, not the deduplicated one.
#
#   bash 05_filter.sh
#
# Input   results/cuttag/bam/<label>.rmDup.bam
# Output  results/cuttag/bam/<label>.filt.bam  (+ .bai)
#
# MIN_MAPQ is set in config/config.sh, default 30. Secondary and supplementary
# alignments are dropped as well, so each fragment is counted once.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

BAM="${OUT_DIR}/cuttag/bam"
MAPQ="${MIN_MAPQ:-30}"
mkdir -p "${BAM}"

COUNT_TSV="${BAM}/filtered_read_counts.tsv"
printf 'label\treads\n' > "${COUNT_TSV}"

found=0
for dt in cuttag atac; do
    while IFS=$'\t' read -r cond rep label mark input; do
        in_bam="${BAM}/${label}.rmDup.bam"
        if [[ ! -f "${in_bam}" ]]; then
            echo "[skip] ${label}: ${in_bam} not found (run 04 first)" >&2
            continue
        fi
        found=$((found + 1))
        out_bam="${BAM}/${label}.filt.bam"

        echo "[${label}] samtools view -q ${MAPQ}" >&2
        samtools view -b -F 0x900 -q "${MAPQ}" -o "${out_bam}" "${in_bam}"
        samtools index "${out_bam}"

        printf '%s\t%s\n' "${label}" \
            "$(samtools view -c -F 0x900 "${out_bam}")" \
            >> "${COUNT_TSV}"
    done < <(samples "${dt}")
done

if (( found == 0 )); then
    echo "ERROR: no deduplicated BAM under ${BAM}. Run 04_remove_dup.sh first." >&2
    exit 1
fi

echo "done -> ${BAM}/"
echo "read counts -> ${COUNT_TSV}"
