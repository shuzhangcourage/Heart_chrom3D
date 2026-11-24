# Heart_chrom3D

Analysis code accompanying the manuscript on **chromatin 3D organisation during
postnatal mouse heart development and in Cdca5-overexpressing cardiomyocytes**.

Four conditions are compared throughout:

| `cond` | Meaning | Other labels used in raw data |
|---|---|---|
| `day1` | postnatal day 1 | `N1`, `P1` |
| `day7` | postnatal day 7 | `N7`, `P7` |
| `day56` | postnatal day 56 | `N56`, `P56` |
| `cdca5oe` | Cdca5 overexpression | `cdca5`, `C227` |

Data types: **Hi-C** (8 libraries, 2 biological replicates per condition),
**CUT&Tag**, **ATAC-seq** and **RNA-seq**.

Genome: `mm10`. CUT&Tag and ATAC-seq libraries were spiked in with *E. coli* and
aligned against both genomes.

---

## Repository layout

| Path | Contents |
|---|---|
| `config/` | `init.sh`, `config.example.sh`, `select.sh`, `samples.tsv`, `contrasts.tsv` |
| `env/` | Conda environment and software versions |
| `01_HiC/` | Hi-C: alignment, matrices, QC, decay, compartments, TADs, boundary classification, loops |
| `02_CUTandTag_ATAC/` | CUT&Tag and ATAC-seq, one shared alignment pipeline distinguished by the `datatype` column |
| `03_RNAseq/` | STAR → featureCounts → DESeq2 |
| `docs/` | Workflow, figure-to-code map, naming map, references |
| `data/` | How to obtain input data; no data files are stored here |

Each directory README lists its analysis steps with the tools and outputs, the
run order, and the parameters that matter.

## Quick start

```bash
git clone https://github.com/shuzhangcourage/Heart_chrom3D.git
cd Heart_chrom3D

# 1. environment
conda env create -f env/environment.yml
conda activate heart_chrom3d

# 2. configuration (config/config.sh is git-ignored)
cp config/config.example.sh config/config.sh
$EDITOR config/config.sh

# 3. point samples.tsv at your own files, then run a step
bash 01_HiC/<step>.sh
```

Every script starts with `source "$(dirname "$0")/../config/init.sh"`, which
loads `config/config.sh` (paths) and `config/select.sh` (sample selection), and
takes no positional arguments. A step is re-run simply by editing
`config/samples.tsv`.

### How conditions are iterated

No script hard-codes a sample name. Scripts ask `config/select.sh` for the
subset they need:

```bash
source config/init.sh

while IFS=$'\t' read -r cond rep label mark input; do
    echo "processing ${cond} ${label}"
done < <(samples hic)          # datatype: hic | cuttag | atac | chipseq | rnaseq
```

Adding a sample means adding one line to `config/samples.tsv`; no script
changes. `config/contrasts.tsv` does the same job for DESeq2 comparisons.

---

## Workflow

[`docs/workflow.md`](docs/workflow.md) gives the full pipeline: which step feeds
which, at what resolution, with which parameters.

## Figure-to-code map

[`docs/figure_legends.md`](docs/figure_legends.md) records, for each panel, the
script, the samples it consumes and the intermediate files it was drawn from.

## Naming

Three naming conventions coexist in the raw data (`day1` / `N1` / `P1` for the
same condition, and `H3K27` → `k27ac`, `H3K9` → `k9me3`). They are reconciled in
exactly one place, `config/samples.tsv`. See
[`docs/naming_map.md`](docs/naming_map.md).

## Data availability

Raw and processed sequencing data: *accession to be added on publication*.
[`data/README.md`](data/README.md) documents the expected directory layout and
the reference files that are not distributed with this repository.

Public ChIP-seq tracks used in the Cdca5 correlation analysis (H3K27me3, H3K4me1,
H3K4me3, H3K9me3) are listed in `config/samples.tsv` with `source = public`.

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite the manuscript once published.

## License

MIT — see [`LICENSE`](LICENSE).

## Contact

Shu Zhang · Guangzhou Medical University · `shuzhangcourage` on GitHub
