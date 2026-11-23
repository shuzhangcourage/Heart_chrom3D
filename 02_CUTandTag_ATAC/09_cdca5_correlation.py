#!/usr/bin/env python3
"""
09_cdca5_correlation.py

Spearman correlation between binned Cdca5 signal and the chromatin marks.

    python3 09_cdca5_correlation.py \
        --matrix results/correlation/counts_100kb.npz \
        --meta   results/correlation/tracks.tsv \
        --outdir results/correlation

multiBigwigSummary writes an .npz with one row of per-bin averages per bigWig.
Correlation is computed on ranks (Spearman), because the distributions of
CUT&Tag, ATAC-seq and ChIP-seq signal are heavy tailed and far from normal;
Pearson on raw coverage would mostly report the few very strong bins.

Bins where any track has no value are dropped before ranking, so all pairs are
compared over the same set of genomic bins.

Two figures are produced:
  correlation_matrix.pdf   every track against every other track
  cdca5_correlation.pdf    per condition, Cdca5 against each mark measured there
"""

import argparse
import csv
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def read_meta(path):
    meta = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            meta[row["label"]] = (row["cond"], row["mark"])
    return meta


def load_matrix(path):
    data = np.load(path, allow_pickle=True)
    matrix = np.asarray(data["matrix"], dtype=float)
    labels = [str(x) for x in data["labels"]]
    if matrix.ndim != 2:
        sys.exit(f"ERROR: unexpected matrix shape {matrix.shape} in {path}")
    if len(labels) != matrix.shape[0]:
        sys.exit(
            f"ERROR: {len(labels)} labels but {matrix.shape[0]} rows in {path}"
        )
    return matrix, labels


def rank_rows(matrix):
    """Column-wise ranks of every row, NaNs already removed by the caller."""
    order = np.argsort(matrix, axis=1, kind="mergesort")
    ranks = np.empty(matrix.shape, dtype=float)
    rows = np.arange(matrix.shape[0])[:, None]
    ranks[rows, order] = np.arange(matrix.shape[1], dtype=float)[None, :]
    return ranks


def spearman(matrix):
    finite = np.isfinite(matrix).all(axis=0)
    kept = matrix[:, finite]
    if kept.shape[1] < 3:
        sys.exit(
            "ERROR: fewer than three bins are shared by all tracks. "
            "Check that the bigWigs cover the same genome."
        )
    ranks = rank_rows(kept)
    # constant tracks would divide by zero; guard them explicitly
    std = ranks.std(axis=1)
    if (std == 0).any():
        bad = np.where(std == 0)[0]
        sys.exit(f"ERROR: track(s) with constant signal: {bad.tolist()}")
    corr = np.corrcoef(ranks)
    print(f"bins used: {kept.shape[1]} of {matrix.shape[1]}", file=sys.stderr)
    return corr


def write_table(corr, labels, path):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["track"] + labels)
        for label, row in zip(labels, corr):
            writer.writerow([label] + [f"{v:.4f}" for v in row])


def plot_matrix(corr, labels, path):
    size = max(6.0, 0.42 * len(labels) + 3.0)
    fig, ax = plt.subplots(figsize=(size, size))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticklabels(labels, fontsize=7)

    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(
                j, i, f"{corr[i, j]:.2f}",
                ha="center", va="center", fontsize=5,
                color="white" if abs(corr[i, j]) > 0.6 else "black",
            )

    ax.set_title("Spearman correlation, 100 kb bins", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def cdca5_pairs(corr, labels, meta):
    """For each condition, Cdca5 against every other mark measured there."""
    by_cond = {}
    for idx, label in enumerate(labels):
        cond, mark = meta.get(label, ("unknown", "unknown"))
        by_cond.setdefault(cond, {})[mark] = idx

    pairs = []
    for cond, marks in by_cond.items():
        if "Cdca5" not in marks:
            continue
        ref = marks["Cdca5"]
        for mark, idx in marks.items():
            if mark == "Cdca5":
                continue
            pairs.append((cond, mark, corr[ref, idx]))
    return pairs


def plot_cdca5(pairs, path):
    if not pairs:
        print(
            "WARNING: no Cdca5 track paired with another mark in the same "
            "condition; skipping cdca5_correlation.pdf",
            file=sys.stderr,
        )
        return

    conds = []
    for cond, _, _ in pairs:
        if cond not in conds:
            conds.append(cond)

    width = 0.8 / max(len(conds), 1)
    fig, ax = plt.subplots(figsize=(max(5.0, 1.1 * len({p[1] for p in pairs}) + 2), 4))
    marks = sorted({p[1] for p in pairs})
    positions = np.arange(len(marks))

    for k, cond in enumerate(conds):
        values = [next((v for c, m, v in pairs if c == cond and m == mark), np.nan)
                  for mark in marks]
        ax.bar(positions + k * width, values, width, label=cond)

    ax.set_xticks(positions + width * (len(conds) - 1) / 2)
    ax.set_xticklabels(marks, rotation=45, ha="right")
    ax.set_ylabel("Spearman rho with Cdca5")
    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_ylim(-1, 1)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Cdca5 correlation with chromatin marks", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Spearman correlation of binned Cdca5 and mark signal."
    )
    parser.add_argument("--matrix", required=True, help="multiBigwigSummary npz")
    parser.add_argument("--meta", required=True, help="label -> condition, mark")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    matrix, labels = load_matrix(args.matrix)
    meta = read_meta(args.meta)
    corr = spearman(matrix)

    write_table(corr, labels, os.path.join(args.outdir, "correlation.tsv"))
    plot_matrix(corr, labels, os.path.join(args.outdir, "correlation_matrix.pdf"))
    plot_cdca5(
        cdca5_pairs(corr, labels, meta),
        os.path.join(args.outdir, "cdca5_correlation.pdf"),
    )
    print(f"done -> {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
