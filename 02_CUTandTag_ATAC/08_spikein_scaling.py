#!/usr/bin/env python3
"""
08_spikein_scaling.py

Turns raw read counts into per-library spike-in scaling factors.

    python3 08_spikein_scaling.py \
        --counts results/cuttag/spikein_counts.tsv \
        --output results/cuttag/spikein_factors.tsv

The idea behind spike-in normalisation: an equal amount of E. coli chromatin is
added to every library before sequencing, so the number of reads that align to
the E. coli genome measures how much material actually went into the tube and
how efficiently it was tagmented. Two libraries with the same biology but
different spike-in counts are therefore not comparable directly.

    scale_factor = median(spike-in reads) / spike-in reads of this library

A library with half the median spike-in signal gets a factor of 2. Multiplying
its CPM coverage track by that factor brings it back onto the same scale as the
others, and the resulting tracks can be compared across conditions.

The median is used rather than the mean so that one failed library cannot drag
every factor with it. Libraries with zero spike-in reads cannot be rescaled and
are reported with an empty factor; 06_bigwig.sh skips them.
"""

import argparse
import csv
import statistics
import sys


def read_counts(path):
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
    return rows


def compute_factors(rows):
    counts = []
    for row in rows:
        try:
            value = int(row["spikein_reads"])
        except (KeyError, TypeError, ValueError):
            value = 0
        row["spikein_reads_int"] = value
        if value > 0:
            counts.append(value)

    if not counts:
        sys.exit(
            "ERROR: no library has any spike-in read. Check that 02_mapping.sh\n"
            "       aligned to the E. coli genome and that SPIKEIN_INDEX is right."
        )

    baseline = statistics.median(counts)
    for row in rows:
        spikein = row["spikein_reads_int"]
        mm10 = int(row.get("mm10_reads") or 0)
        total = mm10 + spikein
        row["spikein_fraction"] = f"{spikein / total:.6f}" if total else "NA"
        row["scale_factor"] = f"{baseline / spikein:.4f}" if spikein else ""

    return rows, baseline


def write_factors(rows, baseline, path):
    fields = [
        "label",
        "scale_factor",
        "cond",
        "rep",
        "mark",
        "mm10_reads",
        "spikein_reads",
        "spikein_fraction",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"spike-in baseline (median reads): {baseline:.0f}", file=sys.stderr)


def print_summary(rows, baseline):
    print(f"# spike-in baseline = {baseline:.0f} reads", file=sys.stderr)
    print(f"{'library':<20}{'spike-in':>12}{'fraction':>10}{'factor':>10}",
          file=sys.stderr)
    for row in rows:
        factor = row["scale_factor"] or "n/a"
        print(
            f"{row['label']:<20}{row['spikein_reads_int']:>12}"
            f"{row['spikein_fraction']:>10}{factor:>10}",
            file=sys.stderr,
        )
    low = [r["label"] for r in rows if not r["scale_factor"]]
    if low:
        print(
            "WARNING: no spike-in reads for: " + ", ".join(low),
            file=sys.stderr,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Compute spike-in scaling factors from read counts."
    )
    parser.add_argument("--counts", required=True, help="spikein_counts.tsv")
    parser.add_argument("--output", required=True, help="spikein_factors.tsv")
    args = parser.parse_args()

    rows = read_counts(args.counts)
    if not rows:
        sys.exit(f"ERROR: {args.counts} is empty.")

    rows, baseline = compute_factors(rows)
    write_factors(rows, baseline, args.output)
    print_summary(rows, baseline)


if __name__ == "__main__":
    main()
