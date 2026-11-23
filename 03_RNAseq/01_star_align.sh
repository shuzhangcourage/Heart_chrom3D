#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 03_RNAseq/01_star_align.sh
#
# STAR alignment for the RNA-seq libraries.
#
#   bash 01_star_align.sh
#
# Input   <R1 fastq> and <R2 fastq>, path of R1 taken from the "input" column
#         of config/samples.tsv for rows with datatype = rnaseq
# Output  results/rnaseq/bam/<label>.bam (+ .bai)
#         results/rnaseq/star/<label>.Log.final.out
#
# The STAR index is generated once into $STAR_INDEX if it is not there yet;
# generating it takes roughly as much memory as the genome, so MEM_MB is read
# from config/config.sh.
#
# R2 is derived from R1 by replacing the mate tag in the file name: _R1 -> _R2,
# _1 -> _2, or R1 -> R2, whichever file actually exists.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

OUT="${OUT_DIR}/rnaseq"
BAM="${OUT}/bam"
LOG="${OUT}/star"
mkdir -p "${BAM}" "${LOG}"

derive_r2() {
    local r1="$1" cand
    for cand in "${r1/_R1/_R2}" "${r1/_1./_2.}" "${r1/R1/R2}"; do
        if [[ "${cand}" != "${r1}" && -f "${cand}" ]]; then
            printf '%s\n' "${cand}"
            return 0
        fi
    done
    return 1
}

# --- index, once ------------------------------------------------------------
if [[ ! -f "${STAR_INDEX}/SAindex" ]]; then
    echo "STAR index not found, generating ${STAR_INDEX}" >&2
    mkdir -p "${STAR_INDEX}"
    STAR --runMode genomeGenerate \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${REF_FA}" \
        --sjdbGTFfile "${GTF}" \
        --sjdbOverhang 100 \
        --limitGenomeGenerateRAM "${MEM_MB}"000000
fi

# --- align ------------------------------------------------------------------
found=0
while IFS=$'\t' read -r cond rep label mark input; do
    r1="${DATA_DIR}/${input}"
    if [[ ! -f "${r1}" ]]; then
        echo "[skip] ${label}: ${r1} not found" >&2
        continue
    fi
    if ! r2="$(derive_r2 "${r1}")"; then
        echo "[skip] ${label}: mate 2 of ${r1} not found" >&2
        continue
    fi
    found=$((found + 1))

    echo "[${label}] STAR" >&2
    STAR --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${r1}" "${r2}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${LOG}/${label}." \
        --outBAMsortingThreadN "${THREADS}"

    mv "${LOG}/${label}.Aligned.sortedByCoord.out.bam" "${BAM}/${label}.bam"
    samtools index "${BAM}/${label}.bam"
done < <(samples rnaseq)

if (( found == 0 )); then
    echo "ERROR: no RNA-seq fastq pair found. Check rows with datatype = rnaseq" >&2
    echo "       in config/samples.tsv: the input column must point at R1." >&2
    exit 1
fi

echo "done -> ${BAM}/"
