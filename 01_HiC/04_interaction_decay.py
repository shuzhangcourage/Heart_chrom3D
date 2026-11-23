#!/usr/bin/env python3
"""
Intra-chromosomal interaction decay (contact probability vs genomic separation).

For each sample the cis expected contact frequency is computed with cooltools,
smoothed, and the local log-log slope is derived. The figure has two panels:
the decay curve and its slope, which makes the position of the plateau and the
shoulder comparable between conditions.

Original code
    scripts/01_Interaction_Decay.py (day1 / day7 / day56 hard-coded at 5 kb)
Changes here
    - samples, resolution, view and output come from the command line
    - the cooler path may be an .mcool or a single-resolution .cool
    - values are written to a TSV alongside the figure
    - the dist < 2 mask and the drop_duplicates('dist') aggregation of the
      original are kept, so the curve matches the published figure

Usage
    python 04_interaction_decay.py \
        --sample day1=/path/to/day1.mcool \
        --sample day7=/path/to/day7.mcool \
        --resolution 5000 --out interaction_decay.pdf
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

import cooler  # noqa: E402
import cooltools  # noqa: E402

PALETTE = ["#1f4e79", "#c0504d", "#4f8b3b", "#8064a2", "#e07b39", "#31859c"]


def open_cooler(path: str, resolution: int) -> cooler.Cooler:
    """Open either an .mcool (at the requested resolution) or a .cool file."""
    p = str(path)
    if p.endswith(".mcool"):
        return cooler.Cooler(f"{p}::/resolutions/{resolution}")
    return cooler.Cooler(p)


def read_view(path: str | None) -> pd.DataFrame | None:
    """Read a BED3+ view file (chrom, start, end[, name])."""
    if not path:
        return None
    df = pd.read_csv(path, sep="\t", header=None)
    df = df.iloc[:, :4]
    if df.shape[1] == 3:
        df[3] = df[0] + ":" + df[1].astype(str) + "-" + df[2].astype(str)
    df.columns = ["chrom", "start", "end", "name"]
    return df


def decay_table(clr, view, smooth_sigma: float, threads: int,
                min_dist_bins: int) -> pd.DataFrame:
    cvd = cooltools.expected_cis(
        clr=clr,
        view_df=view,
        smooth=True,
        aggregate_smoothed=True,
        smooth_sigma=smooth_sigma,
        nproc=threads,
    )
    col = "balanced.avg.smoothed.agg"
    cvd.loc[cvd["dist"] < min_dist_bins, col] = np.nan

    tab = cvd.drop_duplicates(subset=["dist"])[["dist_bp", col]].copy()
    tab = tab.dropna()
    tab = tab.sort_values("dist_bp").reset_index(drop=True)
    tab["slope"] = np.gradient(np.log(tab[col].to_numpy()),
                               np.log(tab["dist_bp"].to_numpy()))
    return tab.rename(columns={col: "contact_frequency"})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", action="append", required=True,
                    metavar="NAME=PATH",
                    help="condition name and cooler path; repeat per sample")
    ap.add_argument("--resolution", type=int, default=5000)
    ap.add_argument("--view", default=None,
                    help="BED file of chromosomes or chromosome arms")
    ap.add_argument("--smooth-sigma", type=float, default=0.1)
    ap.add_argument("--min-dist-bins", type=int, default=2)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--xmin", type=float, default=1e4)
    ap.add_argument("--xmax", type=float, default=1e8)
    ap.add_argument("--out", required=True, help="output figure (pdf/png)")
    ap.add_argument("--out-table", default=None, help="output TSV of the curves")
    args = ap.parse_args()

    view = read_view(args.view)

    tables: dict[str, pd.DataFrame] = {}
    for spec in args.sample:
        if "=" not in spec:
            ap.error(f"--sample expects NAME=PATH, got {spec!r}")
        name, path = spec.split("=", 1)
        print(f"[{name}] {path}", file=sys.stderr)
        clr = open_cooler(path, args.resolution)
        tables[name] = decay_table(clr, view, args.smooth_sigma,
                                   args.threads, args.min_dist_bins)

    if not tables:
        print("no samples could be read", file=sys.stderr)
        return 1

    fig, axs = plt.subplots(
        figsize=(6.5, 9),
        nrows=2,
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )
    for i, (name, tab) in enumerate(tables.items()):
        colour = PALETTE[i % len(PALETTE)]
        axs[0].loglog(tab["dist_bp"], tab["contact_frequency"], "-",
                      color=colour, label=name)
        axs[1].semilogx(tab["dist_bp"], tab["slope"], "-",
                        color=colour, alpha=0.8)

    axs[0].set_ylabel("IC contact frequency")
    axs[0].set_xlim(args.xmin, args.xmax)
    axs[0].grid(False)
    axs[0].legend(frameon=False, loc="upper right")
    axs[1].set_xlabel("separation, bp")
    axs[1].set_ylabel("slope")
    axs[1].grid(False)
    fig.tight_layout()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print(f"figure  -> {out}", file=sys.stderr)

    if args.out_table:
        long = pd.concat(
            [tab.assign(sample=name) for name, tab in tables.items()],
            ignore_index=True,
        )[["sample", "dist_bp", "contact_frequency", "slope"]]
        out_tab = Path(args.out_table)
        out_tab.parent.mkdir(parents=True, exist_ok=True)
        long.to_csv(out_tab, sep="\t", index=False)
        print(f"table   -> {out_tab}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
