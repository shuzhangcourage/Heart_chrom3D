#!/usr/bin/env python3
"""
Library statistics table (Table 1).

Reads pairtools statistics files and writes a table with, per library: total
read pairs, unmapped, mapped, PCR duplicates, duplicates removed, cis, trans,
valid pairs, and the cis distance fractions when they are available.

Where the numbers come from
    The counts are read from the pairtools stats files written by step 01
    (<label>.parse.stats.txt and <label>.dedup.stats.txt). Rates are computed
    here. The cis distance fractions are only filled in if the stats files were
    produced with pairtools distance filters, which adds keys such as
    "cis_1kb+" and "cis_10kb+"; otherwise those cells stay empty.

Usage
    python 03_qc_table1.py \
        --sample N1_rep1=/path/N1_rep1.dedup.stats.txt \
        --out results/hic/qc/table1.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

COLUMNS = [
    "library",
    "total_pairs",
    "unmapped",
    "mapped",
    "mapped_rate",
    "pcr_duplicates",
    "pcr_duplicate_rate",
    "duplicates_removed",
    "duplicates_removed_rate",
    "cis",
    "cis_rate",
    "trans",
    "trans_rate",
    "valid_pairs",
    "valid_rate",
    "cis_lt_1kb",
    "cis_ge_1kb",
    "cis_ge_10kb",
]


def read_stats(path: str) -> dict[str, float]:
    """Parse a pairtools stats file: one 'key<TAB>value' pair per line."""
    stats: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                stats[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return stats


def rate(num: float | None, den: float | None) -> float | None:
    if num is None or den in (None, 0):
        return None
    return round(num / den, 4)


def build_row(library: str, stats: dict[str, float]) -> dict:
    total = stats.get("total")
    unmapped = stats.get("total_unmapped")
    mapped = stats.get("total_mapped")
    dups = stats.get("total_dups")
    nodups = stats.get("total_nodups")
    cis = stats.get("cis")
    trans = stats.get("trans")

    row = {
        "library": library,
        "total_pairs": total,
        "unmapped": unmapped,
        "mapped": mapped,
        "mapped_rate": rate(mapped, total),
        "pcr_duplicates": dups,
        "pcr_duplicate_rate": rate(dups, mapped),
        "duplicates_removed": nodups,
        "duplicates_removed_rate": rate(nodups, total),
        "cis": cis,
        "cis_rate": rate(cis, nodups),
        "trans": trans,
        "trans_rate": rate(trans, nodups),
        "valid_pairs": (cis + trans) if (cis is not None and trans is not None) else None,
        "valid_rate": rate((cis + trans) if (cis is not None and trans is not None)
                           else None, total),
        "cis_lt_1kb": stats.get("cis_<1kb"),
        "cis_ge_1kb": stats.get("cis_1kb+"),
        "cis_ge_10kb": stats.get("cis_10kb+"),
    }
    return row


def to_int_columns(df: pd.DataFrame) -> pd.DataFrame:
    counts = ["total_pairs", "unmapped", "mapped", "pcr_duplicates",
              "duplicates_removed", "cis", "trans", "valid_pairs",
              "cis_lt_1kb", "cis_ge_1kb", "cis_ge_10kb"]
    for c in counts:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce").round().astype("Int64")
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", action="append", required=True,
                    metavar="LABEL=STATS_FILE")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--out-markdown", default=None,
                    help="also write a markdown version of the table")
    args = ap.parse_args()

    rows = []
    for spec in args.sample:
        if "=" not in spec:
            ap.error("--sample expects LABEL=STATS_FILE")
        label, path = spec.split("=", 1)
        if not Path(path).exists():
            print(f"[skip] {label}: {path} not found", file=sys.stderr)
            continue
        rows.append(build_row(label, read_stats(path)))

    if not rows:
        print("no statistics files could be read", file=sys.stderr)
        return 1

    df = pd.DataFrame(rows).reindex(columns=COLUMNS)
    df = to_int_columns(df)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"table -> {out}", file=sys.stderr)

    if args.out_markdown:
        md = Path(args.out_markdown)
        md.parent.mkdir(parents=True, exist_ok=True)
        md.write_text(df.to_markdown(index=False) + "\n")
        print(f"markdown -> {md}", file=sys.stderr)

    print(df.to_string(index=False), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
