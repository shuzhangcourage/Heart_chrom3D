#!/usr/bin/env python3
"""
Inter-chromosome contact heatmap (trans expected contacts).

The average balanced contact frequency between every pair of chromosomes is
pivoted into a square matrix and drawn as a log-scaled heatmap. Because each
chromosome is compared with every other, this is also called the trans
interaction heatmap.

Original code
    scripts/02_trans_heatmap.py (day1 only, 5 kb, chromosome list hard-coded)
Changes here
    - the sample, resolution, chromosome list and output come from the CLI
    - the matrix is written to a TSV so the figure can be redrawn or restyled
    - the title uses the sample name instead of the literal "Day1"

Note on colours
    The matrix is divided by its maximum and then drawn with LogNorm, as in the
    original. That keeps the appearance of the published figure; switch to
    --no-normalise if absolute frequencies are wanted.

Usage
    python 05_interchrom_heatmap.py \
        --sample day1=/path/to/day1.mcool --resolution 5000 \
        --out interchrom_day1.pdf
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
from matplotlib import colors  # noqa: E402

import cooler  # noqa: E402
import cooltools  # noqa: E402

DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 20)]


def open_cooler(path: str, resolution: int) -> cooler.Cooler:
    """Open either an .mcool (at the requested resolution) or a .cool file."""
    p = str(path)
    if p.endswith(".mcool"):
        return cooler.Cooler(f"{p}::/resolutions/{resolution}")
    return cooler.Cooler(p)


def trans_matrix(clr, chroms: list[str], threads: int) -> pd.DataFrame:
    # keep only chromosomes that are actually present
    present = [c for c in chroms if c in set(clr.chromnames)]
    missing = [c for c in chroms if c not in present]
    if missing:
        print(f"warning: chromosomes not in the cooler, skipped: "
              f"{','.join(missing)}", file=sys.stderr)
    if len(present) < 2:
        raise SystemExit("need at least two chromosomes for a trans matrix")

    view = pd.DataFrame({
        "chrom": present,
        "start": 0,
        "end": [clr.chromsizes[c] for c in present],
        "name": present,
    })

    ac = cooltools.expected_trans(clr, view_df=view, nproc=threads)
    mat = (ac.pivot_table(values="balanced.avg",
                          index="region1",
                          columns="region2",
                          observed=True)
             .reindex(index=present, columns=present))
    return mat


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True, metavar="NAME=PATH")
    ap.add_argument("--resolution", type=int, default=5000)
    ap.add_argument("--chroms", default=None,
                    help="comma separated list, default chr1..chr19")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cmap", default="YlOrRd")
    ap.add_argument("--no-normalise", action="store_true",
                    help="do not divide the matrix by its maximum")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-table", default=None)
    args = ap.parse_args()

    if "=" not in args.sample:
        ap.error("--sample expects NAME=PATH")
    name, path = args.sample.split("=", 1)

    chroms = ([c.strip() for c in args.chroms.split(",")]
              if args.chroms else DEFAULT_CHROMS)

    clr = open_cooler(path, args.resolution)
    mat = trans_matrix(clr, chroms, args.threads)

    labels = list(mat.index)
    values = mat.to_numpy(dtype=float)

    if not args.no_normalise:
        vmax = np.nanmax(values)
        if vmax > 0:
            values = values / vmax

    # the matrix should be symmetric; fill any gaps from the transpose
    values = np.where(np.isnan(values), values.T, values)

    width = max(6.0, 0.42 * len(labels))
    fig, axs = plt.subplots(
        figsize=(width, width * 0.92),
        ncols=2,
        gridspec_kw={"width_ratios": [20, 1], "wspace": 0.08},
    )
    ax, cax = axs

    im = ax.imshow(values, cmap=args.cmap, norm=colors.LogNorm(), aspect=1.0)
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=10)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_title(args.title or f"{name} average interactions between chromosomes",
                 fontsize=12)

    ax.set_xticks([x - 0.5 for x in range(1, len(labels))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(labels))], minor=True)
    ax.grid(which="minor", linewidth=0.4, color="0.75")

    fig.colorbar(im, cax=cax,
                 label="IC contact frequency"
                       + ("" if args.no_normalise else " (scaled to max)"))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print(f"figure -> {out}", file=sys.stderr)

    if args.out_table:
        mat.to_csv(args.out_table, sep="\t")
        print(f"table  -> {args.out_table}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
