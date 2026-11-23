#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 02_CUTandTag_ATAC/08_spikein_scaling.sh
#
# Counts spike-in reads and turns them into per-library scaling factors.
#
#   bash 08_spikein_scaling.sh
#
# Input   results/cuttag/sam/<label>_spikein.sam   aligned to E. coli
#         results/cuttag/bam/<label>.filt.bam      aligned to mm10
# Output  results/cuttag/spikein_counts.tsv         raw counts
#         results/cuttag/spikein_factors.tsv        factor per library
#
# The factor is computed by 08_spikein_scaling.py. Feed the factors back into
# the coverage step to get tracks that are comparable across libraries:
#
#   SPIKEIN_TABLE=results/cuttag/spikein_factors.tsv bash 06_bigwig.sh
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

SAM="${OUT_DIR}/cuttag/sam"
BAM="${OUT_DIR}/cuttag/bam"
RAW="${OUT_DIR}/cuttag/spikein_counts.tsv"
MAPQ="${MIN_MAPQ:-30}"

printf 'cond\trep\tlabel\tmark\tmm10_reads\tspikein_reads\n' > "${RAW}"

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    sp_sam="${SAM}/${label}_spikein.sam"
    in_bam="${BAM}/${label}.filt.bam"
    if [[ ! -f "${sp_sam}" ]]; then
        echo "[skip] ${label}: ${sp_sam} not found (run 02 first)" >&2
        continue
    fi
    if [[ ! -f "${in_bam}" ]]; then
        echo "[skip] ${label}: ${in_bam} not found (run 05 first)" >&2
        continue
    fi
    found=$((found + 1))

    echo "[${label}] counting reads" >&2
    mm10_reads="$(samtools view -c -F 0x900 "${in_bam}")"
    spikein_reads="$(samtools view -c -F 0x900 -q "${MAPQ}" "${sp_sam}")"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${cond}" "${rep}" "${label}" "${mark}" \
        "${mm10_reads}" "${spikein_reads}" >> "${RAW}"
done < <(samples cuttag)

if (( found == 0 )); then
    echo "ERROR: no spike-in alignment found under ${SAM}." >&2
    echo "       Run 02_mapping.sh first." >&2
    exit 1
fi

python3 "$(dirname "$0")/08_spikein_scaling.py" \
    --counts "${RAW}" \
    --output "${OUT_DIR}/cuttag/spikein_factors.tsv"

echo "done -> ${OUT_DIR}/cuttag/spikein_factors.tsv"
