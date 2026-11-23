#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 03_RNAseq/03_deseq2.R
#
# Differential expression between conditions.
#
#   Rscript 03_deseq2.R
#
# Input   results/rnaseq/counts/counts_matrix.tsv   from step 02
#         config/samples.tsv                        condition of each library
#         config/contrasts.tsv                      which comparisons to run
# Output  results/rnaseq/deseq2/normalised_counts.tsv
#         results/rnaseq/deseq2/<contrast>.tsv       one table per comparison
#         results/rnaseq/deseq2/summary.tsv         up/down counts per comparison
#
# The design is ~condition. Comparisons are not hard-coded: every row of
# config/contrasts.tsv becomes one results table, so comparisons can be added
# or removed without touching this script.
#
# Log2 fold changes are shrunk with apeglm when the package is installed and
# the coefficient exists; otherwise the unshrunken results are written.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(DESeq2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
repo_dir <- normalizePath(file.path(dirname(script_path), ".."))

out_dir <- Sys.getenv("OUT_DIR", file.path(repo_dir, "results"))
counts_file <- file.path(out_dir, "rnaseq", "counts", "counts_matrix.tsv")
samples_file <- file.path(repo_dir, "config", "samples.tsv")
contrasts_file <- file.path(repo_dir, "config", "contrasts.tsv")
res_dir <- file.path(out_dir, "rnaseq", "deseq2")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# --- inputs -----------------------------------------------------------------
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
samples <- read.delim(samples_file, comment.char = "#", stringsAsFactors = FALSE)
samples <- samples[samples$datatype == "rnaseq", , drop = FALSE]

missing <- setdiff(colnames(counts), samples$label)
if (length(missing) > 0) {
    stop("libraries in the count matrix but not in samples.tsv: ",
         paste(missing, collapse = ", "))
}
samples <- samples[match(colnames(counts), samples$label), , drop = FALSE]

if (nrow(samples) == 0) {
    stop("no rows with datatype = rnaseq in ", samples_file)
}

cat("libraries:", nrow(samples), "\n")
print(table(samples$cond))

col_data <- data.frame(
    condition = factor(samples$cond),
    row.names = samples$label
)

# --- model ------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = col_data,
    design = ~condition
)

# genes with almost no reads only slow DESeq2 down and add noise
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("genes kept:", nrow(dds), "of", nrow(counts), "\n")

dds <- DESeq(dds)

norm <- counts(dds, normalized = TRUE)
write.table(
    data.frame(gene_id = rownames(norm), norm, check.names = FALSE),
    file.path(res_dir, "normalised_counts.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

# --- comparisons ------------------------------------------------------------
contrasts <- read.delim(contrasts_file, comment.char = "#", stringsAsFactors = FALSE)

summary_rows <- list()
for (i in seq_len(nrow(contrasts))) {
    name <- contrasts$contrast[i]
    numerator <- contrasts$numerator[i]
    denominator <- contrasts$denominator[i]

    coef_name <- paste0("condition_", numerator, "_vs_", denominator)
    if (requireNamespace("apeglm", quietly = TRUE) &&
        coef_name %in% resultsNames(dds)) {
        res <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
    } else {
        res <- results(dds, contrast = c("condition", numerator, denominator))
    }

    tbl <- data.frame(
        gene_id = rownames(res),
        baseMean = res$baseMean,
        log2FoldChange = res$log2FoldChange,
        lfcSE = res$lfcSE,
        stat = res$stat,
        pvalue = res$pvalue,
        padj = res$padj,
        check.names = FALSE
    )
    tbl <- tbl[order(tbl$padj, -abs(tbl$log2FoldChange)), ]
    write.table(
        tbl,
        file.path(res_dir, paste0(name, ".tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
    )

    sig <- !is.na(tbl$padj) & tbl$padj < 0.05 & abs(tbl$log2FoldChange) > 1
    summary_rows[[name]] <- data.frame(
        contrast = name,
        tested = sum(!is.na(tbl$pvalue)),
        significant = sum(sig),
        up = sum(sig & tbl$log2FoldChange > 0),
        down = sum(sig & tbl$log2FoldChange < 0),
        stringsAsFactors = FALSE
    )
    cat(sprintf("%-22s up %5d  down %5d\n", name,
                summary_rows[[name]]$up, summary_rows[[name]]$down))
}

summary_tbl <- do.call(rbind, summary_rows)
write.table(
    summary_tbl,
    file.path(res_dir, "summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

cat("done ->", res_dir, "\n")
