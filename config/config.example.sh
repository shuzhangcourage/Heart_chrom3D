#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# config/config.example.sh
#
# Copy to config/config.sh and edit. config/config.sh is git-ignored, so local
# paths never end up in the repository.
#
#   cp config/config.example.sh config/config.sh
#
# Every script in this repository starts with
#   source "$(dirname "$0")/../config/config.sh"
# so all paths below become available as environment variables.
# ---------------------------------------------------------------------------

# --- repository root --------------------------------------------------------
export PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- input data -------------------------------------------------------------
# Where the published data were downloaded to. See data/README.md.
export DATA_DIR="${PROJECT_DIR}/data"

export HIC_DIR="${DATA_DIR}/hic"          # *.mcool / *.cool
export BAM_DIR="${DATA_DIR}/bam"          # CUT&Tag, ATAC-seq, RNA-seq BAM
export BW_DIR="${DATA_DIR}/bigwig"        # coverage tracks *.bw
export PEAK_DIR="${DATA_DIR}/peaks"       # MACS2 output
export PUBLIC_DIR="${DATA_DIR}/public"    # public ChIP-seq tracks

# --- reference genome (mm10) ------------------------------------------------
export REF_DIR="${DATA_DIR}/reference"
export GENOME="${REF_DIR}/mm10.chrom.sizes"
export REF_FA="${REF_DIR}/mm10.fa"
export GTF="${REF_DIR}/mm10.gtf"
export BLACKLIST="${REF_DIR}/mm10_blacklist_v2.bed"

# STAR index directory (RNA-seq)
export STAR_INDEX="${REF_DIR}/star_mm10"

# E. coli spike-in reference (CUT&Tag / ATAC-seq)
export SPIKEIN_INDEX="${REF_DIR}/ecoli/ecoli"

# --- output -----------------------------------------------------------------
export OUT_DIR="${PROJECT_DIR}/results"

# --- compute ----------------------------------------------------------------
export THREADS=16
export MEM_MB=32000

# --- tool specific defaults -------------------------------------------------
# optimalTAD only uses --chipseq to pick the best window size; h3k4me3 or
# k27ac are both acceptable. The final window size is fixed at 40000.
export OPTIMALTAD_CHIP="${PUBLIC_DIR}/h3k4me3.bw"

# mustache loop calling
export MUSTACHE_PT=0.01

# MACS2 uses a q-value cutoff
export MACS2_QVALUE=0.05

# min mapping quality kept for CUT&Tag / ATAC-seq BAM
export MIN_MAPQ=30

mkdir -p "${OUT_DIR}"
