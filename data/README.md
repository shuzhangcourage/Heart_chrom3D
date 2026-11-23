# Data

**No data files are stored in this repository.** This directory only documents
how to obtain and lay out the inputs.

## Expected layout

Point `DATA_DIR` in `config/config.sh` at a directory with this structure. All
paths in `config/samples.tsv` are relative to it.

```
data/
├── hic/          *.mcool (ICE-balanced, 1 kb base resolution)
├── bam/          CUT&Tag / ATAC-seq / RNA-seq BAM after filtering
├── bigwig/       coverage tracks (*.bw)
├── peaks/        MACS2 narrowPeak output
├── public/       public ChIP-seq tracks used for comparison
└── reference/    mm10.fa, mm10.gtf, mm10.chrom.sizes, mm10_blacklist_v2.bed,
                  star_mm10/, ecoli/ecoli.*
```

## Processed data (this study)

*Deposited on submission — GEO accession pending.*

When the accession is available, replace this paragraph and update the
"Data availability" section of the top-level README.

## Reference genome

`mm10`. Exact source and release are listed in `env/versions.md`.
The E. coli spike-in genome supplied with the CUT&Tag kit is used as a second
alignment reference for CUT&Tag and ATAC-seq libraries.

## Public data

The correlation analysis (`02_CUTandTag_ATAC/09_cdca5_correlation.py`) compares
Cdca5 CUT&Tag signal with H3K27me3, H3K4me1, H3K4me3 and H3K9me3 tracks.
Rows with `source = public` in `config/samples.tsv` point to these tracks.

Source and accession: *to be completed on publication.*

## Files that are too large for Git

All of the above. `.gitignore` excludes `data/` and every common genomic file
extension; only this README is tracked.
