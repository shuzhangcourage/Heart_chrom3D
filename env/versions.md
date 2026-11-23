# Software versions

Versions recorded for the analyses in this repository. Tools not listed are
given in the manuscript methods; `env/environment.yml` pins the versions used to
build the conda environment.

## Hi-C 3.0

| Tool | Version |
|---|---|
| BWA | 0.7.19-r1273 |
| pairtools | 1.1.3 |
| cooler | 0.10.2 |
| cooltools | 0.7.1 |
| Juicer Tools | 1.22.01 (`.hic` generation only) |
| optimalTAD | 0.1.0 |
| coolpup.py / plotpup.py | 1.1.0 |
| mustache | 1.3.3 |

## CUT&Tag, ATAC-seq, RNA-seq

| Tool | Version |
|---|---|
| bowtie2 | *see methods* |
| picard | *see methods* |
| samtools | *see methods* |
| deepTools | *see methods* |
| MACS2 | *see methods* |
| STAR | *see methods* |
| featureCounts (subread) | *see methods* |
| DESeq2 | *see methods* |

## Reference files

| File | Source |
|---|---|
| `mm10.fa` | UCSC mm10 |
| `mm10.gtf` | Ensembl, release used for the alignment |
| `mm10.chrom.sizes` | derived from `mm10.fa` |
| `mm10_blacklist_v2.bed` | ENCODE mm10 blacklist |
| E. coli spike-in genome | reference supplied with the CUT&Tag kit |

Download locations are documented in `data/README.md`.

## Parameter notes

**optimalTAD.** `--resolution 10000 --window_size_min 40000
--window_size_max 40000`. The `--chipseq` track only selects the optimal window
size; H3K4me3 and H3K27ac give the same result, and the track in use is set by
`$OPTIMALTAD_CHIP`.

**mustache.** `-r 5kb -pt 0.01`, restricted to `chr1`–`chr19`.

**Inter-chromosome heatmap.** Computed at 5 kb; the same script accepts
`--resolution 500000` for a 500 kb rendering.

**MACS2.** q-value cutoff, `$MACS2_QVALUE` (default 0.05).
