#!/usr/bin/env python3
"""
Boundary refinement: turn TAD domain edges into one boundary bin each.

optimalTAD reports a boundary as a single coordinate. This step picks the bin
that actually carries the boundary:

  1. for every TAD, take the insulation score at its start and at its end
  2. drop bins without a score
  3. sort by chromosome and position
  4. when two neighbouring bins are only one bin apart, keep the one with the
     lower insulation score and drop the other
  5. keep autosomes chr1-chr19

Original code
    script 01-extract_all_anchor_insulation.py
    script 02-boudary_insulation_score_filtered.py
    script 03-autosome_insulation_score.py
    merged here; the day1-only assumption is gone

Note
    The original also dropped boundaries with a positive insulation score. That
    filter is exposed as --drop-positive (off by default) because the published
    boundary set was also edited by hand afterwards; enabling it reproduces that
    intermediate state but not the hand edits.

Usage
    python 11_boundary_refine_IS.py \
        --sample day1=/path/day1_IS.tsv:/path/day1.domains.tad \
        --bin-size 10000 --out-dir results/hic/boundary
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def read_is(path: str, window: int) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    score = f"log2_insulation_score_{window}"
    if score not in df.columns:
        candidates = [c for c in df.columns if "insulation_score" in c]
        if not candidates:
            raise SystemExit(f"{path}: no insulation score column found")
        score = candidates[0]
        print(f"  using column {score}", file=sys.stderr)
    return df[["chrom", "start", score]].rename(columns={score: "score"})


def read_domains(path: str) -> pd.DataFrame:
    # optimalTAD writes .domains.tad comma separated; accept tabs too
    df = pd.read_csv(path, sep=None, engine="python")
    cols = {c.lower(): c for c in df.columns}
    need = [cols.get("chr") or cols.get("chrom"), cols.get("start"), cols.get("end")]
    if any(c is None for c in need):
        raise SystemExit(f"{path}: expected columns Chr, Start, End")
    df = df.rename(columns={need[0]: "chrom", need[1]: "start", need[2]: "end"})
    return df[["chrom", "start", "end"]]


def refine(is_df: pd.DataFrame, domains: pd.DataFrame,
           bin_size: int, drop_positive: bool) -> pd.DataFrame:
    index = is_df.set_index(["chrom", "start"])["score"]

    starts = domains[["chrom", "start"]].rename(columns={"start": "pos"})
    ends = domains[["chrom", "end"]].rename(columns={"end": "pos"})
    anchors = pd.concat([starts, ends], ignore_index=True)

    def lookup(row):
        try:
            return index.loc[(row["chrom"], row["pos"])]
        except KeyError:
            return float("nan")

    anchors["score"] = anchors.apply(lookup, axis=1)
    anchors = anchors.rename(columns={"pos": "start"})
    anchors = anchors.dropna(subset=["score"])
    anchors = anchors.sort_values(["chrom", "start"]).reset_index(drop=True)

    kept = []
    i = 0
    starts_arr = anchors["start"].to_numpy()
    chroms = anchors["chrom"].to_numpy()
    scores = anchors["score"].to_numpy()
    while i < len(anchors):
        if (i < len(anchors) - 1
                and chroms[i] == chroms[i + 1]
                and abs(starts_arr[i + 1] - starts_arr[i]) == bin_size):
            kept.append(i if scores[i] <= scores[i + 1] else i + 1)
            i += 2
        else:
            kept.append(i)
            i += 1

    out = anchors.iloc[kept].copy()
    if drop_positive:
        out = out[out["score"] <= 0]

    out = out[out["chrom"].astype(str).str.match(r"^chr([1-9]|1[0-9])$")]
    out["end"] = out["start"] + bin_size
    return out[["chrom", "start", "end", "score"]].reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", action="append", required=True,
                    metavar="LABEL=IS_TSV:DOMAINS_TAD")
    ap.add_argument("--bin-size", type=int, default=10000)
    ap.add_argument("--window", type=int, default=40000,
                    help="insulation window used in step 10")
    ap.add_argument("--drop-positive", action="store_true")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for spec in args.sample:
        if "=" not in spec:
            ap.error("--sample expects LABEL=IS_TSV:DOMAINS_TAD")
        label, rest = spec.split("=", 1)
        parts = [p for p in rest.split(":") if p]
        if len(parts) != 2:
            ap.error(f"--sample {label}: need IS_TSV:DOMAINS_TAD")
        is_path, tad_path = parts

        if not Path(is_path).exists() or not Path(tad_path).exists():
            print(f"[skip] {label}: input missing", file=sys.stderr)
            continue

        is_df = read_is(is_path, args.window)
        domains = read_domains(tad_path)
        bed = refine(is_df, domains, args.bin_size, args.drop_positive)

        tsv = out_dir / f"{label}_boundary.tsv"
        bed[["chrom", "start", "end", "score"]].to_csv(tsv, sep="\t", index=False)
        bed3 = out_dir / f"{label}_boundary.bed"
        bed[["chrom", "start", "end"]].to_csv(bed3, sep="\t", header=False, index=False)

        print(f"[{label}] {len(domains)} TADs -> {len(bed)} boundaries "
              f"({tsv.name})", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
