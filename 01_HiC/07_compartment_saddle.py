#!/usr/bin/env python3
"""
Compartment A/B saddle plot and compartment strength.

The saddle plot is the average observed/expected contact frequency between bins
grouped by their eigenvector value: a corner-to-corner pattern means the two
compartments are well separated, a flat plot means they are not.

Original code
    scripts/06-cdca5_plot.py (single sample, values hard-coded)
Changes here
    - samples, inputs and output paths come from the command line
    - the plotting function is kept as in the original, only guarded against a
      missing qrange
    - compartment strength is written to a table and drawn as a bar chart

Strength
    The trimmed saddle matrix is split at the bin where the binned eigenvector
    changes sign. With A = eigenvector > 0 and B = eigenvector < 0,

        strength = (mean(AA) + mean(BB)) / (mean(AB) + mean(BA))

    which is the same quantity as the AA+BB over AB+BA ratio computed from a
    matrix dump, but read off the saddle matrix directly.

Usage
    python 07_compartment_saddle.py \
        --sample day1=/path/day1.mcool,/path/day1_eigs.tsv,/path/day1_expected.tsv \
        --resolution 100000 --out saddle.pdf
"""

from __future__ import annotations

import argparse
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import ticker  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec  # noqa: E402

import cooler  # noqa: E402
import cooltools  # noqa: E402


# ---------------------------------------------------------------------------
# plotting, adapted from the original script
# ---------------------------------------------------------------------------
def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = {vmin % 10 * 10, vmax % 10, 1}

        def __call__(self, x, pos=None):
            if x not in (vmin, 1, vmax):
                return ""
            return "{x:g}".format(x=x)

    if qrange is None:
        qrange = (0.0, 1.0)
    lo, hi = qrange

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    groupmean = track[track.columns[3]].groupby(
        digitized_track[digitized_track.columns[3]]
    ).mean()

    binedges = np.linspace(lo, hi, n_bins + 1)

    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    if subplot_spec is not None:
        gs_factory = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    else:
        gs_factory = GridSpec

    grid = {}
    gs = gs_factory(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = {**fig_kws_default, **(fig_kws or {})}
        fig = plt.figure(**fig_kws)

    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("only linear and log colour scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws = {**dict(cmap=cmap, rasterized=True), **(heatmap_kws or {})}
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    margin_kws = {**dict(edgecolor="k", facecolor=color, linewidth=1),
                  **(margin_kws or {})}

    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])
    plt.barh(binedges, height=1 / len(binedges), width=groupmean,
             align="edge", **margin_kws)
    plt.xlim(plt.xlim()[1], plt.xlim()[0])
    plt.ylim(hi, lo)
    for side in ("top", "bottom", "left"):
        plt.gca().spines[side].set_visible(False)
    plt.gca().xaxis.set_visible(False)

    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])
    plt.bar(binedges, width=1 / len(binedges), height=groupmean,
            align="edge", **margin_kws)
    plt.xlim(lo, hi)
    for side in ("top", "right", "left"):
        plt.gca().spines[side].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws = {**dict(fraction=0.8, label=clabel or ""), **(cbar_kws or {})}
    if scale == "linear" and vmin is not None and vmax is not None:
        cb = plt.colorbar(img, **cbar_kws)
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        cb = plt.colorbar(img, format=MinOneMaxFormatter(),
                          cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid["ax_heatmap"].grid(False)
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    return grid


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------
def open_cooler(path: str, resolution: int) -> cooler.Cooler:
    if path.endswith(".mcool"):
        return cooler.Cooler(f"{path}::/resolutions/{resolution}")
    return cooler.Cooler(path)


def autosome_view(clr) -> pd.DataFrame:
    """chr1..chr19, or the bare names if the cooler has no 'chr' prefix."""
    import re

    chroms = [c for c in clr.chromnames if re.fullmatch(r"chr([1-9]|1[0-9])", c)]
    if not chroms:
        chroms = [c for c in clr.chromnames if re.fullmatch(r"([1-9]|1[0-9])", c)]
    sizes = clr.chromsizes.loc[chroms]
    return pd.DataFrame({"chrom": chroms, "start": 0,
                         "end": sizes.values, "name": chroms})


def load_eigenvector(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    e1 = cols.get("e1")
    if e1 is None:
        raise SystemExit(f"{path}: no E1 column")
    return df[["chrom", "start", "end", e1]].rename(columns={e1: "E1"})


def compartment_strength(saddledata: np.ndarray, n_bins: int,
                         qrange: tuple[float, float]) -> tuple[float, float, float]:
    """AA+BB over AB+BA on the trimmed saddle matrix."""
    C = saddledata
    if (C.shape[0] - n_bins) == 2:
        C = C[1:-1, 1:-1]
    binedges = np.linspace(qrange[0], qrange[1], n_bins + 1)
    split = int(np.clip(np.searchsorted(binedges, 0.0), 1, n_bins - 1))

    with np.errstate(invalid="ignore"):
        aa = np.nanmean(C[split:, split:])
        bb = np.nanmean(C[:split, :split])
        ab = np.nanmean(C[:split, split:])
        ba = np.nanmean(C[split:, :split])

    intra = aa + bb
    inter = ab + ba
    strength = float(intra / inter) if inter else float("nan")
    return strength, float(intra), float(inter)


def parse_sample(spec: str) -> tuple[str, str, str, str]:
    if "=" not in spec:
        raise SystemExit("--sample expects NAME=COOL,EIGS,EXPECTED")
    name, rest = spec.split("=", 1)
    parts = [p for p in rest.split(",") if p]
    if len(parts) != 3:
        raise SystemExit(
            f"--sample {name}: need COOLER,EIGS,EXPECTED, got {len(parts)} field(s)")
    return name, parts[0], parts[1], parts[2]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", action="append", required=True,
                    metavar="NAME=COOL,EIGS,EXPECTED")
    ap.add_argument("--resolution", type=int, default=100000)
    ap.add_argument("--n-bins", type=int, default=50)
    ap.add_argument("--q-lo", type=float, default=0.025)
    ap.add_argument("--q-hi", type=float, default=0.975)
    ap.add_argument("--vmin", type=float, default=0.5)
    ap.add_argument("--vmax", type=float, default=2.0)
    ap.add_argument("--scale", choices=["log", "linear"], default="log")
    ap.add_argument("--cmap", default="coolwarm")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True, help="single-sample figure prefix")
    ap.add_argument("--out-stacked", default=None,
                    help="figure with all samples side by side")
    ap.add_argument("--out-strength", default=None,
                    help="TSV with the compartment strength per sample")
    args = ap.parse_args()

    qrange = (args.q_lo, args.q_hi)
    cache = {}
    rows = []

    for spec in args.sample:
        name, cool_path, eig_path, exp_path = parse_sample(spec)
        print(f"[{name}] {cool_path}", file=sys.stderr)

        clr = open_cooler(cool_path, args.resolution)
        eigvec = load_eigenvector(eig_path)
        cvd = pd.read_csv(exp_path, sep="\t")

        total, count = cooltools.saddle(
            clr, cvd, eigvec, "cis",
            n_bins=args.n_bins, qrange=qrange, view_df=autosome_view(clr),
        )
        saddledata = total / count

        strength, intra, inter = compartment_strength(saddledata, args.n_bins,
                                                      qrange)
        rows.append({"sample": name, "strength": strength,
                     "intra": intra, "inter": inter})

        cache[name] = (eigvec, saddledata)

        plt.figure(figsize=(5, 5))
        saddleplot(
            eigvec, saddledata, args.n_bins,
            qrange=qrange, vmin=args.vmin, vmax=args.vmax,
            cmap=args.cmap, scale=args.scale,
            clabel="average observed/expected contact frequency",
            title=f"{name} {args.resolution // 1000}kb cis saddle (chr1-19)",
            xlabel="E1", ylabel="E1",
        )
        stem = f"{args.out}_{name}" if len(args.sample) > 1 else args.out
        Path(stem).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(f"{stem}.pdf", bbox_inches="tight")
        plt.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
        plt.close("all")
        print(f"  saddle  -> {stem}.pdf", file=sys.stderr)

    if args.out_stacked and len(cache) > 1:
        fig = plt.figure(figsize=(5 * len(cache), 5.5),
                         constrained_layout=True)
        parent = fig.add_gridspec(nrows=1, ncols=len(cache))
        for i, (name, (eigvec, saddledata)) in enumerate(cache.items()):
            saddleplot(
                eigvec, saddledata, args.n_bins,
                qrange=qrange, vmin=args.vmin, vmax=args.vmax,
                cmap=args.cmap, scale=args.scale,
                clabel=None, title=name, subplot_spec=parent[0, i],
            )
        fig.suptitle(f"cis saddle, chr1-19, {args.resolution // 1000}kb",
                     fontsize=14, y=1.02)
        Path(args.out_stacked).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out_stacked, bbox_inches="tight")
        print(f"stacked -> {args.out_stacked}", file=sys.stderr)
        plt.close("all")

    table = pd.DataFrame(rows)
    print(table.to_string(index=False), file=sys.stderr)

    if args.out_strength:
        Path(args.out_strength).parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.out_strength, sep="\t", index=False)
        print(f"strength-> {args.out_strength}", file=sys.stderr)

        fig, ax = plt.subplots(figsize=(max(4, 0.9 * len(table)), 3.6))
        ax.bar(table["sample"], table["strength"],
               color="#1f4e79", width=0.6)
        ax.set_ylabel("compartment strength\n(AA+BB)/(AB+BA)")
        ax.set_ylim(0, max(1.0, float(np.nanmax(table["strength"])) * 1.25))
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="x", rotation=30)
        fig.tight_layout()
        bar = str(Path(args.out_strength).with_suffix(".pdf"))
        fig.savefig(bar, bbox_inches="tight")
        print(f"bar     -> {bar}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
