#!/usr/bin/env python3
"""Plot saved truth-side INDEL recall by length and repeat overlap."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[1]
LENGTHS = ("1-5", "6-15", "16+")
CALLERS = ("deepvariant", "gatk")
COLORS = {"deepvariant": "#0072B2", "gatk": "#D55E00"}
LABELS = {"deepvariant": "DeepVariant 1.0.0", "gatk": "GATK 4.1.8.0"}
PANELS = (
    ("repeat_overlap", "Overlapping tandem repeats"),
    ("no_repeat_overlap", "Not overlapping tandem repeats"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", default="MYT1")
    parser.add_argument("--input-dir")
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    gene = args.gene.upper()

    if not gene or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in gene):
        raise ValueError("Invalid gene name")

    source = (
        Path(args.input_dir).resolve() / "summary.csv"
        if args.input_dir
        else ROOT / f"results/context_analysis/{gene}/summary.csv"
    )
    rows = {}

    with source.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["caller"], row["length_group"], row["repeat_context"])
            if key in rows:
                raise ValueError(f"Duplicate summary row: {key}")
            if row["gene"] != gene:
                raise ValueError("Gene label does not match requested gene")

            total = int(row["truth_total"])
            tp, fn = int(row["truth_tp"]), int(row["truth_fn"])
            if min(total, tp, fn) < 0 or tp + fn != total:
                raise ValueError(f"Invalid counts: {key}")

            if total:
                if abs(float(row["recall"]) - tp / total) > 1e-9:
                    raise ValueError(f"Recall does not match counts: {key}")
            elif row["recall"] != "":
                raise ValueError(f"Empty group has nonempty recall: {key}")

            rows[key] = (total, tp)

    expected = {
        (caller, length, context)
        for caller in CALLERS
        for length in LENGTHS
        for context, _ in PANELS
    }
    if set(rows) != expected:
        raise ValueError("Missing or unexpected summary groups")

    for length in LENGTHS:
        for context, _ in PANELS:
            if rows[("deepvariant", length, context)][0] != rows[("gatk", length, context)][0]:
                raise ValueError("Truth denominators differ between callers")

    output = Path(args.output_dir).resolve() if args.output_dir else ROOT / "results/figures"
    targets = [
        output / f"{gene}_indel_recall_by_context.{ext}"
        for ext in ("png", "svg", "pdf")
    ]
    if any(path.exists() for path in targets):
        raise FileExistsError("Figure already exists; nothing overwritten")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
    })

    fig, axes = plt.subplots(1, 2, figsize=(11, 6), sharey=True)
    width = 0.32

    for ax, (context, title) in zip(axes, PANELS):
        ax.set_title(title, fontsize=11, pad=14)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E3E7EC", linewidth=0.8)
        ax.set_ylim(0, 112)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.yaxis.set_major_formatter(PercentFormatter(100))
        ax.set_xticks(range(3), ["1–5", "6–15", "16+"])
        ax.set_xlabel("Maximum called INDEL length per truth record (bp)")
        ax.set_xlim(-0.6, 2.6)

        for i, length in enumerate(LENGTHS):
            total = rows[("deepvariant", length, context)][0]
            if total == 0:
                ax.text(i, 5, "No records\n(n = 0)",
                        ha="center", va="bottom", color="#5B6470", fontsize=9)
                continue

            for j, caller in enumerate(CALLERS):
                total, tp = rows[(caller, length, context)]
                recall = 100 * tp / total
                x = i + (j - 0.5) * width

                ax.bar(x, recall, width=width,
                       color=COLORS[caller], zorder=3)
                ax.text(x, recall + 2, f"{tp}/{total}",
                        ha="center", va="bottom", fontsize=9)

    axes[0].set_ylabel("INDEL recall")
    fig.suptitle(
        f"{gene}: INDEL recall by length and tandem-repeat overlap",
        fontsize=15, fontweight="bold", y=0.97
    )
    fig.legend(
        handles=[Patch(facecolor=COLORS[c], label=LABELS[c]) for c in CALLERS],
        loc="upper center", bbox_to_anchor=(0.5, 0.915),
        ncol=2, frameon=False
    )
    fig.text(
        0.07, 0.055,
        "HG002 · GRCh38 · GIAB v5.0q evaluable regions. Labels show TP / truth records.\n"
        "Repeat overlap uses the REF span (insertion anchor). Length uses the largest called literal allele change.\n"
        "Exploratory record-level analysis; observations may be linked within haplotype blocks.",
        fontsize=8.5, color="#46515F", va="bottom"
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.77,
                        bottom=0.25, wspace=0.16)

    output.mkdir(parents=True, exist_ok=True)
    for path in targets:
        fig.savefig(path, dpi=300, facecolor="white")
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
