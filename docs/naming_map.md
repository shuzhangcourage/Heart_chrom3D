# Naming map

Three naming conventions coexist in the raw data. They are reconciled in
exactly one place — `config/samples.tsv` — and nowhere else.

## Conditions

| Meaning | Hi-C libraries | CUT&Tag / ATAC libraries | `cond` used in this repo |
|---|---|---|---|
| Postnatal day 1 | `N1_rep1`, `N1_rep2` | `P1_*` | `day1` |
| Postnatal day 7 | `N7_rep1`, `N7_rep2` | `P7_*` | `day7` |
| Postnatal day 56 | `N56_rep1`, `N56_rep2` | `P56_*` | `day56` |
| Cdca5 overexpression, clone C227 | `cdca5_rep1`, `cdca5_rep2` | `C227_*` | `cdca5oe` |

`cdca5`, `C227` and `cdca5oe` all refer to the same condition. The clone number
C227 is recorded in the `label` column, not in `cond`.

## Marks

| Raw label | `mark` in this repo | Note |
|---|---|---|
| `CTCF` | `CTCF` | |
| `RAD21`, `Rad21` | `RAD21` | capitalisation varies between files |
| `H3K27` | `k27ac` | H3K27ac |
| `H3K9` | `k9me3` | H3K9me3 |
| `ATAC` | `ATAC` | |
| `Cdca5` | `Cdca5` | Cdca5 CUT&Tag |
| `k27me3`, `k4me1`, `k4me3` | unchanged | public ChIP-seq tracks |

## Genome-wide rules applied during curation

1. Every `cond` value is one of `day1`, `day7`, `day56`, `cdca5oe`.
2. Every `mark` value uses the spelling in the table above.
3. Filenames appear only in `label` and `input`; they are never parsed in a
   script to recover a condition or a mark.
4. Replicate-merged libraries carry `rep = merged`.
