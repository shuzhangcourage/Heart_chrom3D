# 03_RNAseq

STAR → featureCounts → DESeq2.

## Analysis steps

| # | Step | Tool | Output |
|---|---|---|---|
| 01 | alignment | STAR | coordinate-sorted BAM |
| 02 | gene-level read counting | featureCounts | count matrix |
| 03 | differential expression | DESeq2 | results per contrast |

## Run order

`01 → 02 → 03`. All samples are aligned and counted into a single matrix; the
comparisons are then defined once in `config/contrasts.tsv`.

## Parameters

**01 STAR**

```
# index, once
STAR --runMode genomeGenerate --runThreadN $THREADS \
     --genomeDir $STAR_INDEX --genomeFastaFiles $REF_FA \
     --sjdbGTFfile $GTF --sjdbOverhang 100

# align
STAR --runThreadN $THREADS --genomeDir $STAR_INDEX \
     --readFilesIn <R1.fq.gz> <R2.fq.gz> --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix $OUT_DIR/<label>.
```

**02 featureCounts**

```
featureCounts -T $THREADS -p --countReadPairs \
              -t exon -g gene_id -a $GTF \
              -o counts_matrix.txt <bam files...>
```

`-p --countReadPairs` counts fragments rather than reads and assumes paired-end
libraries.

**03 DESeq2**

`~condition`, with the comparisons listed in `config/contrasts.tsv`. The default
is all six pairwise contrasts between `day1`, `day7`, `day56` and `cdca5oe`.

```
dds  <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition)
dds  <- DESeq(dds)
res  <- results(dds, contrast = c("condition", numerator, denominator))
```

## Notes

- The design formula can be extended with a batch term if the libraries were
  prepared in more than one batch.
- Log2 fold changes are shrunk with apeglm before plotting.
- Comparisons come from `config/contrasts.tsv`; delete rows that are not needed.
