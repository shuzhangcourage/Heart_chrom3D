#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_HiC/01_mapping_pairtools.sh
#
# Align Hi-C 3.0 reads and produce pairs files.
#
# Input naming convention
#   ${FASTQ_DIR}/hic/<label>_R1.fq.gz and _R2.fq.gz
#   FASTQ_DIR defaults to ${DATA_DIR}/fastq
#
# Output under results/hic/pairs/
#   <label>.pairs             deduplicated pairs, sorted
#   <label>.parse.stats.txt   pairtools parse statistics
#   <label>.dedup.stats.txt   pairtools dedup statistics
#   <label>.PT.bam / .bai     pairs as a name-sorted BAM, for inspection
#
# Parameters follow the manuscript: bwa mem -5SP -T0, min MAPQ 40, walks-policy
# 5unique, max inter-align gap 30.
#
# Duplicates
#   --mark-dups is kept because this is the command that produced the published
#   matrices. Note what it does: duplicate pairs are flagged rather than
#   dropped, so they are still in the .pairs file that reaches cooler and they
#   are counted in the matrix. Remove the flag to get a strictly deduplicated
#   matrix; the duplicate counts are written to the stats file either way.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

FASTQ_DIR="${FASTQ_DIR:-${DATA_DIR}/fastq}"
OUT="${OUT_DIR}/hic/pairs"
TMPDIR_="${TMPDIR:-/tmp}"
mkdir -p "${OUT}" "${TMPDIR_}/pairtools"

found=0
while IFS=$'\t' read -r cond rep label mark input; do
    r1="${FASTQ_DIR}/hic/${label}_R1.fq.gz"
    r2="${FASTQ_DIR}/hic/${label}_R2.fq.gz"
    if [[ ! -f "${r1}" || ! -f "${r2}" ]]; then
        echo "[skip] ${label}: ${r1} or ${r2} not found" >&2
        continue
    fi
    found=$((found + 1))
    echo "[${label}] aligning" >&2

    bwa mem -5SP -T0 -t "${THREADS}" "${REF_FA}" "${r1}" "${r2}" \
    | pairtools parse \
        --min-mapq 40 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --nproc-in "${THREADS}" \
        --nproc-out "${THREADS}" \
        --chroms-path "${GENOME}" \
        --output-stats "${OUT}/${label}.parse.stats.txt" \
    | pairtools sort \
        --tmpdir "${TMPDIR_}/pairtools" \
        --nproc "${THREADS}" \
    | pairtools dedup \
        --nproc-in "${THREADS}" \
        --nproc-out "${THREADS}" \
        --mark-dups \
        --output-stats "${OUT}/${label}.dedup.stats.txt" \
    | pairtools split \
        --nproc-in "${THREADS}" \
        --nproc-out "${THREADS}" \
        --output-pairs "${OUT}/${label}.pairs" \
        --output-sam - \
    | samtools view -bS -@ "${THREADS}" - \
    | samtools sort -@ "${THREADS}" -o "${OUT}/${label}.PT.bam" -

    samtools index "${OUT}/${label}.PT.bam"
done < <(samples hic)

if (( found == 0 )); then
    echo "ERROR: no Hi-C fastq found under ${FASTQ_DIR}/hic" >&2
    echo "       expected <label>_R1.fq.gz / _R2.fq.gz for each label" >&2
    echo "       in config/samples.tsv" >&2
    exit 1
fi

echo "done -> ${OUT}/"
