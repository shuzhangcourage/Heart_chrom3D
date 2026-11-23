# Workflow

Every step reads `config/samples.tsv` through `config/select.sh`; no condition
or filename is hard-coded. Steps are independent shell scripts that can be run
on their own once their inputs exist.

---

## 1. Hi-C 3.0 (`01_HiC/`)

| # | Script | In | Out |
|---|---|---|---|
| 01 | `01_mapping_pairtools.sh` | paired fastq | `*.pairs` |
| 02 | `02_matrix_cooler.sh` | `*.pairs` | `*_matrix_1kb.mcool` |
| 03 | `03_qc_table1.py` | pairtools stats, mcool | Table 1 |
| 04 | `04_interaction_decay.py` | mcool | intra-chromosomal decay curve |
| 05 | `05_interchrom_heatmap.py` | mcool | inter-chromosome contact heatmap |
| 06 | `06_compartment_eigs.sh` | mcool | E1 eigenvectors (50 kb) |
| 07 | `07_compartment_saddle.py` | mcool + E1 | saddle plot, compartment strength |
| 08 | `08_compartment_strength.py` | saddle output | strength per condition |
| 09 | `09_compartment_switch.py` | E1 tracks | A/B switch positions, day1 vs day56 |
| 10 | `10_optimalTAD.sh` | mcool + ChIP track | TAD domains |
| 11 | `11_boundary_refine_IS.py` | domains + insulation score | refined boundary BED |
| 12 | `12_boundary_IS_track.sh` | mcool | IS bigWig, 1st and 2nd derivative bigWigs |
| 13 | `13_boundary_matrix.sh` | IS derivative bigWig + boundaries | matrix ±200 kb, relative values, k-means |
| 14 | `14_boundary_IS_heatmap.sh` | matrix | IS heatmaps per cluster |
| 15 | `15_boundary_clusters_pileup.sh` | cluster BED + mcool | coolpup pile-ups (±500 kb) |
| 16 | `16_loops_mustache.sh` | mcool | loop calls |

### Alignment parameters (from the methods)

```
bwa mem -5SP -T0
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30
pairtools dedup          # PCR duplicates removed
cooler cload pairix      # + cooler zoomify
```
Resolution: 5 kb for decay, 10 kb for TAD/boundary and pile-ups, 50 kb for
compartments, 100 kb for saddle plots.

### Boundary classification

Follows the approach of *Insulation between adjacent TADs is controlled by the
width of their boundaries* (https://github.com/JSegueni/Boundary_classification_width):

1. insulation score at 40 kb window, 10 kb bins
2. first and second derivative of the IS profile locate the boundary centre
3. matrix of ±200 kb around each boundary, values made relative to the centre bin
4. k-means (k = 3) on the relative matrix, then heatmaps of the real IS per cluster

---

## 2. CUT&Tag and ATAC-seq (`02_CUTandTag_ATAC/`)

Both data types use the same seven alignment steps, following the
[CUT&Tag tutorial](https://yezhengstat.github.io/CUTTag_tutorial/). They differ
only in the `datatype` column of `config/samples.tsv`.

| # | Script | In | Out |
|---|---|---|---|
| 01 | `01_trim_galore.sh` | paired fastq | trimmed fastq |
| 02 | `02_mapping.sh` | trimmed fastq | SAM, mm10 + E. coli spike-in |
| 03 | `03_sort_addRG.sh` | SAM | coordinate-sorted SAM with read group |
| 04 | `04_remove_dup.sh` | sorted SAM | deduplicated BAM |
| 05 | `05_filter.sh` | BAM | MAPQ ≥ 30 BAM |
| 06 | `06_bigwig.sh` | BAM | CPM coverage bigWig |
| 07 | `07_peaks_macs2.sh` | BAM | MACS2 peaks (q-value cutoff) |
| 08 | `08_spikein_scaling.py` | E. coli mapping stats | normalisation factors |
| 09 | `09_cdca5_correlation.py` | bigWigs, 4 conditions × 9 marks, 100 kb bins | Spearman correlation heatmap |

### Alignment parameters (from the methods)

```
trim_galore --paired --retain_unpaired --quality 20 --gzip
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700
picard SortSam SORT_ORDER=coordinate
samtools addreplacerg
picard MarkDuplicates REMOVE_DUPLICATES=true
samtools view -q 30
bamCoverage --binSize 10 --normalizeUsing CPM
```

Peaks: MACS2 with a q-value cutoff (`$MACS2_QVALUE`, default 0.05).

---

## 3. RNA-seq (`03_RNAseq/`)

| # | Script | In | Out |
|---|---|---|---|
| 01 | `01_star_align.sh` | paired fastq | coordinate-sorted BAM |
| 02 | `02_featurecounts.sh` | BAM | gene-level count matrix |
| 03 | `03_deseq2.R` | count matrix | DESeq2 results per contrast |

Comparisons are defined in `config/contrasts.tsv` (all six pairwise contrasts
between the four conditions by default).

---

## 4. Figures

The mapping from figure panel to script is kept in
[`figure_legends.md`](figure_legends.md).
