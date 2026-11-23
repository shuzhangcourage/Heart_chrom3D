# 02_CUTandTag_ATAC

CUT&Tag and ATAC-seq share one alignment pipeline. The commands are identical;
the two data types differ only in the `datatype` column of `config/samples.tsv`.

Alignment follows the [CUT&Tag tutorial](https://yezhengstat.github.io/CUTTag_tutorial/).
Libraries were spiked in with *E. coli*, so every library is aligned to mm10 and
to the spike-in genome separately.

## Analysis steps

| # | Step | Applies to | Output |
|---|---|---|---|
| 01 | adapter trimming | both | trimmed fastq |
| 02 | alignment to mm10 and to the E. coli spike-in | both | SAM |
| 03 | coordinate sorting, read group | both | sorted SAM |
| 04 | duplicate removal | both | deduplicated BAM |
| 05 | mapping-quality filter | both | MAPQ ≥ 30 BAM |
| 06 | coverage track, CPM normalised | both | bigWig |
| 07 | peak calling with a q-value cutoff | both | narrowPeak |
| 08 | spike-in scaling factors | CUT&Tag | normalisation table |
| 09 | Cdca5 correlation with chromatin marks | both + public ChIP-seq | heatmap |

## Run order

`01 → 06` is a single chain per library.
`07` needs the filtered BAM from `05`.
`08` needs the spike-in BAM written by `02`.
`09` needs bigWigs from `06` for all four conditions and all marks, plus the
public ChIP-seq tracks listed in `config/samples.tsv` with `source = public`.

## Parameters

```
trim_galore --paired --retain_unpaired --quality 20 --gzip
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700
picard SortSam SORT_ORDER=coordinate
samtools addreplacerg -r '@RG\tID:RG1\tSM:<sample>\tPL:Illumina\tLB:<library>'
picard MarkDuplicates REMOVE_DUPLICATES=true
samtools view -q 30
bamCoverage --binSize 10 --normalizeUsing CPM
```

Peaks: MACS2 with a q-value cutoff, `$MACS2_QVALUE` (default 0.05),
`-f BAMPE -g mm`.

## Note on step 09

`09` is the only step that reads more than one data type: it bins Cdca5 CUT&Tag
together with CTCF, RAD21, H3K27ac, H3K9me3, H3K27me3, H3K4me1, H3K4me3 and ATAC
into 100 kb bins and computes Spearman correlation per condition. The public
histone tracks come from `config/samples.tsv` with `source = public`.

Mark names are normalised in `config/samples.tsv`: raw `H3K27` becomes `k27ac`
and raw `H3K9` becomes `k9me3`. See [`docs/naming_map.md`](../docs/naming_map.md).
