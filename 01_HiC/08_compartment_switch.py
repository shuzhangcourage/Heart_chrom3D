#!/usr/bin/env python3
"""
Compartment A/B switching between two conditions.

Every bin is labelled A (eigenvector > 0) or B (eigenvector < 0). Bins whose
label differs between the two conditions are reported as switches, split into
B-to-A and A-to-B, and drawn as a genome-wide strip plot.

Orientation
    A positive eigenvector does not by itself mean compartment A: eigs-cis
    returns the eigenvector with an arbitrary sign. Two ways to fix it:
      1. phase it during step 06 with a GC track (--phasing-track); then
         E1 > 0 is A and no signal track is needed here;
      2. pass a signal track per sample (--sample NAME=EIGS,SIGNAL). The sign is
         then flipped if the summed signal is lower on the E1 > 0 side, which is
         the rule used by the earlier A/B assignment script.

Input formats
    EIGS    TSV with chrom, start, end, E1 (the output of step 06)
    SIGNAL  BED/TSV, signal in the last column, one row per bin

Usage
    python 08_compartment_switch.py \
        --sample day1=/path/day1_eigs.tsv \
        --sample day56=/path/day56_eigs.tsv \
        --compare day1=day56 --out switch.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

A_COLOUR = "#c0504d"
B_COLOUR = "#1f4e79"


def load_eigs(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    e1 = cols.get("e1")
    if e1 is None:
        raise SystemExit(f"{path}: no E1 column")
    out = df[["chrom", "start", "end", e1]].rename(columns={e1: "E1"})
    out = out.dropna(subset=["E1"]).reset_index(drop=True)
    return out


def load_signal(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] < 4:
        raise SystemExit(f"{path}: expected at least 4 columns")
    out = df.iloc[:, [0, 1, df.shape[1] - 1]].copy()
    out.columns = ["chrom", "start", "signal"]
    out["signal"] = pd.to_numeric(out["signal"], errors="coerce")
    return out


def orient(eigs: pd.DataFrame, signal: pd.DataFrame | None) -> tuple[pd.DataFrame, bool]:
    """Return the eigenvector with A positive. The flag says whether it flipped."""
    if signal is None:
        return eigs, False

    merged = eigs.merge(signal, on=["chrom", "start"], how="left")
    pos = merged.loc[merged["E1"] > 0, "signal"].sum()
    neg = merged.loc[merged["E1"] < 0, "signal"].sum()
    if not np.isfinite(pos) or not np.isfinite(neg) or pos == neg:
        return eigs, False
    if neg > pos:
        eigs = eigs.copy()
        eigs["E1"] = -eigs["E1"]
        return eigs, True
    return eigs, False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", action="append", required=True,
                    metavar="NAME=EIGS[,SIGNAL]")
    ap.add_argument("--compare", default=None, metavar="LEFT=RIGHT",
                    help="conditions to compare; default: the two given samples")
    ap.add_argument("--min-abs-e1", type=float, default=0.0,
                    help="ignore bins with |E1| below this value")
    ap.add_argument("--chroms", default=None,
                    help="comma separated list, default chr1..chr19")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-bed", default=None, help="switches as a BED file")
    args = ap.parse_args()

    tracks: dict[str, pd.DataFrame] = {}
    for spec in args.sample:
        if "=" not in spec:
            ap.error("--sample expects NAME=EIGS[,SIGNAL]")
        name, rest = spec.split("=", 1)
        parts = [p for p in rest.split(",") if p]
        eigs = load_eigs(parts[0])
        sig = load_signal(parts[1]) if len(parts) > 1 else None
        eigs, flipped = orient(eigs, sig)
        if flipped:
            print(f"[{name}] eigenvector flipped so that A is positive",
                  file=sys.stderr)
        tracks[name] = eigs

    if args.compare:
        if "=" not in args.compare:
            ap.error("--compare expects LEFT=RIGHT")
        left, right = args.compare.split("=", 1)
    elif len(tracks) == 2:
        left, right = list(tracks)[:2]
    else:
        ap.error("give --compare LEFT=RIGHT, or exactly two --sample entries")

    for n in (left, right):
        if n not in tracks:
            ap.error(f"--compare refers to unknown sample {n!r}")

    a = tracks[left][["chrom", "start", "end", "E1"]].rename(columns={"E1": "E1_left"})
    b = tracks[right][["chrom", "start", "end", "E1"]].rename(columns={"E1": "E1_right"})
    m = a.merge(b, on=["chrom", "start", "end"], how="inner")

    if args.min_abs_e1 > 0:
        m = m[(m["E1_left"].abs() >= args.min_abs_e1)
              & (m["E1_right"].abs() >= args.min_abs_e1)]

    m["label_left"] = np.where(m["E1_left"] > 0, "A", "B")
    m["label_right"] = np.where(m["E1_right"] > 0, "A", "B")
    m["switch"] = np.where(
        m["label_left"] != m["label_right"],
        m["label_left"] + "2" + m["label_right"], "-")

    switches = m[m["switch"] != "-"].copy()
    n_ba = int((switches["switch"] == "B2A").sum())
    n_ab = int((switches["switch"] == "A2B").sum())
    print(f"{left} vs {right}: {len(m)} bins compared, "
          f"{len(switches)} switches (B->A {n_ba}, A->B {n_ab})", file=sys.stderr)

    if args.out_bed:
        out_bed = Path(args.out_bed)
        out_bed.parent.mkdir(parents=True, exist_ok=True)
        switches[["chrom", "start", "end", "switch"]].to_csv(
            out_bed, sep="\t", header=False, index=False)
        print(f"bed     -> {out_bed}", file=sys.stderr)

    chroms = ([c.strip() for c in args.chroms.split(",")]
              if args.chroms else [f"chr{i}" for i in range(1, 20)])
    present = [c for c in chroms if c in set(m["chrom"])] or \
              list(dict.fromkeys(m["chrom"]))

    fig, ax = plt.subplots(figsize=(11, 0.32 * len(present) + 1.8))
    for i, chrom in enumerate(present):
        row = m[m["chrom"] == chrom]
        if row.empty:
            continue
        ax.hlines(i, row["start"].min(), row["end"].max(),
                  color="0.9", linewidth=4, zorder=1)
        sw = row[row["switch"] != "-"]
        for label, colour in (("B2A", A_COLOUR), ("A2B", B_COLOUR)):
            sel = sw[sw["switch"] == label]
            if not sel.empty:
                ax.vlines(sel["start"], i - 0.28, i + 0.28,
                          color=colour, linewidth=0.7, zorder=2)

    ax.set_yticks(range(len(present)))
    ax.set_yticklabels(present, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("position, bp")
    ax.set_title(f"Compartment switches, {left} vs {right}\n"
                 f"B->A {n_ba} bins, A->B {n_ab} bins", fontsize=11)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.legend(handles=[Patch(color=A_COLOUR, label="B to A"),
                       Patch(color=B_COLOUR, label="A to B")],
              frameon=False, loc="upper right", ncol=2)
    fig.tight_layout()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print(f"figure  -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
