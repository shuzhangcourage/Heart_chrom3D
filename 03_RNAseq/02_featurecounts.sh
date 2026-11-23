#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 03_RNAseq/02_featurecounts.sh
#
# Gene-level counting for every RNA-seq library into one matrix.
#
#   bash 02_featurecounts.sh
#
# Input   results/rnaseq/bam/<label>.bam
# Output  results/rnaseq/counts/raw_counts.txt        featureCounts output
#         results/rnaseq/counts/counts_matrix.tsv     cleaned, DESeq2 input
#
# -p --countReadPairs counts fragments instead of reads, which is what a
# paired-end library contributes to a gene; -t exon -g gene_id aggregates
# exons into genes. The cleaned matrix keeps gene ids plus one column per
# library in the order of config/samples.tsv.
# ---------------------------------------------------------------------------
set -euo pipefail
source "$(dirname "$0")/../config/init.sh"

BAM="${OUT_DIR}/rnaseq/bam"
CNT="${OUT_DIR}/rnaseq/counts"
mkdir -p "${CNT}"

BAMS=()
LABELS=""
found=0
while IFS=$'\t' read -r cond rep label mark input; do
    bam="${BAM}/${label}.bam"
    if [[ ! -f "${bam}" ]]; then
        echo "[skip] ${label}: ${bam} not found (run 01 first)" >&2
        continue
    fi
    found=$((found + 1))
    BAMS+=("${bam}")
    LABELS="${LABELS:+${LABELS} }${label}"
done < <(samples rnaseq)

if (( found == 0 )); then
    echo "ERROR: no RNA-seq BAM under ${BAM}. Run 01_star_align.sh first." >&2
    exit 1
fi

echo "counting ${found} libraries" >&2
featureCounts \
    -T "${THREADS}" \
    -p --countReadPairs \
    -t exon -g gene_id \
    -a "${GTF}" \
    -o "${CNT}/raw_counts.txt" \
    "${BAMS[@]}"

# keep Geneid plus the count columns, and replace the path header with labels
awk -v OFS='\t' -v labs="${LABELS}" '
    /^#/ { next }
    $1 == "Geneid" {
        printf "Geneid"
        n = split(labs, L, " ")
        for (i = 1; i <= n; i++) printf "%s%s", OFS, L[i]
        printf "\n"
        next
    }
    {
        printf "%s", $1
        for (i = 7; i <= NF; i++) printf "%s%s", OFS, $i
        printf "\n"
    }
' "${CNT}/raw_counts.txt" > "${CNT}/counts_matrix.tsv"

echo "done -> ${CNT}/counts_matrix.tsv"
