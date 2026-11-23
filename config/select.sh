#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# config/select.sh
#
# Sample selection helper. All analysis scripts source this file and query it,
# so that no condition name is ever hard-coded in a script.
#
#   source config/select.sh
#   samples hic                     # all Hi-C records
#   samples cuttag CTCF             # CUT&Tag records for one mark
#   conditions hic                  # distinct conditions for a datatype
#   marks cuttag                    # distinct marks for a datatype
#
# Output is tab separated: cond, rep, label, mark, input
# ---------------------------------------------------------------------------

SAMPLES_TSV="${PROJECT_DIR}/config/samples.tsv"

_read_tsv() {
    # strip comment lines and the column-name header, whatever line they are on
    awk -F'\t' 'NF>1 && !/^#/ && $1!="cond" && $1!="contrast"' "$1"
}

samples() {
    local datatype="$1"
    local mark="${2:-}"
    if [[ -n "${mark}" ]]; then
        _read_tsv "${SAMPLES_TSV}" | awk -F'\t' -v d="${datatype}" -v m="${mark}" \
            '$3==d && $5==m {print $1"\t"$2"\t"$4"\t"$5"\t"$7}'
    else
        _read_tsv "${SAMPLES_TSV}" | awk -F'\t' -v d="${datatype}" \
            '$3==d {print $1"\t"$2"\t"$4"\t"$5"\t"$7}'
    fi
}

conditions() {
    local datatype="$1"
    _read_tsv "${SAMPLES_TSV}" | awk -F'\t' -v d="${datatype}" \
        '$3==d {print $1}' | awk '!seen[$0]++'
}

marks() {
    local datatype="$1"
    _read_tsv "${SAMPLES_TSV}" | awk -F'\t' -v d="${datatype}" \
        '$3==d && $5!="-" {print $5}' | awk '!seen[$0]++'
}

# contrasts for DESeq2: numerator, denominator
contrasts() {
    _read_tsv "${PROJECT_DIR}/config/contrasts.tsv" | awk -F'\t' '{print $1"\t"$2"\t"$3}'
}

export -f samples conditions marks contrasts _read_tsv
