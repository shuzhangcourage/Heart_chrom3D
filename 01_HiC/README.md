# 01_HiC

Hi-C analysis: alignment, matrices, QC, interaction decay, compartments,
TADs, boundary classification and loops.

Every step reads samples from `config/samples.tsv` through `config/select.sh`;
no condition or filename is hard-coded.

## Analysis steps

| # | Step | Tool | Output |
|---|---|---|---|
| 01 | align, parse pairs | bwa, pairtools | `*.pairs` |
| 02 | build matrices, ICE balance, zoomify | cooler | `*_matrix_1kb.mcool` |
| 03 | library statistics | pairtools stats, cooler | Table 1 |
| 04 | intra-chromosomal interaction decay | cooltools | decay curve |
| 05 | inter-chromosome contact heatmap | cooltools | trans heatmap |
| 06 | compartment eigenvectors (100 kb) | cooler, cooltools | E1 tracks |
| 07 | saddle plot and compartment strength (100 kb) | cooltools | saddle plot |
| 08 | compartment switching, day1 vs day56 | bioframe | switch positions |
| 09 | TAD calling (10 kb) | optimalTAD | domains |
| 10 | boundary refinement by insulation score | cooltools | boundary BED |
| 11 | insulation score and its derivatives | deeptools, bedGraphToBigWig | IS + derivative bigWigs |
| 12 | matrix ±200 kb around boundaries, k-means (k = 3) | deeptools | clustered boundaries |
| 13 | insulation score heatmaps per cluster | deeptools | heatmaps |
| 14 | pile-ups per boundary class (±500 kb) | coolpuppy | pile-up plots |
| 15 | loop calling (5 kb) | mustache | loops |

## Run order

`01 → 03` produce the matrices and the QC table.

`04 → 08` (decay and compartments) only need 02.
`09 → 14` form one chain: TAD calling → boundary refinement → insulation score →
derivatives → matrix → clustering → heatmaps → pile-ups.
`15` needs only 02.

The two chains are independent of each other and can run in parallel once the
matrices exist.

## Resolutions

| Analysis | Resolution |
|---|---|
| interaction decay, inter-chromosome heatmap | 5 kb |
| TAD calling, boundary matrix, pile-ups | 10 kb |
| compartments (eigs-cis) | 100 kb |
| saddle plot | 100 kb |

## Parameters worth knowing

**Alignment.** `bwa mem -5SP -T0`; `pairtools parse --min-mapq 40
--walks-policy 5unique --max-inter-align-gap 30`; PCR duplicates removed with
`pairtools dedup`; matrices loaded with `cooler cload pairix` and balanced with
`cooler zoomify --balance`.

**TAD calling.** `--resolution 10000 --window_size_min 40000
--window_size_max 40000`. The `--chipseq` track only selects the optimal window
size and can be either H3K4me3 or H3K27ac; it is exposed as `$OPTIMALTAD_CHIP`.

**Boundary classification.** Follows
[Boundary_classification_width](https://github.com/JSegueni/Boundary_classification_width):
the first and second derivative of the insulation score locate the boundary
centre, a matrix of ±200 kb is made relative to the centre bin, and k-means
(k = 3) separates boundaries by width.

**Relative matrices.** The centre column index is computed from `-b`, `-a` and
the bin size rather than hard-coded, so the matrix extent can be changed safely.

**Pile-ups.** `--flank 500000` at 10 kb; the expected values come from the
`expected_cis` table of the same sample.
