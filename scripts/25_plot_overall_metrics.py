"""Plot overall metrics from the curated benchmark summary tables."""
from pathlib import Path
import csv
import hashlib
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "reports/figures"
PROVENANCE = ROOT / "reports/provenance"
CALLERS = (
    ("deepvariant", "DeepVariant 1.0.0", "#0072B2"),
    ("gatk", "GATK 4.1.8.0", "#D55E00"),
)
METRICS = (
    ("METRIC.Recall", "Recall"),
    ("METRIC.Precision", "Precision"),
    ("METRIC.F1_Score", "F1 score"),
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records = {}
    inputs = []
    for caller, label, color in CALLERS:
        path = ROOT / f"reports/tables/{caller}.summary.csv"
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        inputs.append({
            "path": path.relative_to(ROOT).as_posix(),
            "sha256": sha256(path),
        })
        for variant_type in ("SNP", "INDEL"):
            selected = [
                row for row in rows
                if row["Type"] == variant_type and row["Filter"] == "PASS"
            ]
            if len(selected) != 1:
                raise SystemExit(
                    f"Expected one PASS row: {caller}, {variant_type}"
                )
            row = selected[0]
            for field, _ in METRICS:
                value = float(row[field])
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise SystemExit(f"Invalid metric: {caller}, {field}")
            records[caller, variant_type] = row

    for variant_type in ("SNP", "INDEL"):
        totals = {
            int(records[caller, variant_type]["TRUTH.TOTAL"])
            for caller, _, _ in CALLERS
        }
        if len(totals) != 1 or next(iter(totals)) <= 0:
            raise SystemExit("Truth counts differ or are not positive.")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
    fig.subplots_adjust(
        left=0.09, right=0.98, bottom=0.26, top=0.73, wspace=0.35
    )
    handles = []
    for ax, variant_type in zip(axes, ("SNP", "INDEL")):
        total = int(records["deepvariant", variant_type]["TRUTH.TOTAL"])
        ax.set_title(
            f"{variant_type} | {total} truth records",
            loc="left", fontweight="bold", pad=14,
        )
        for index, (caller, label, color) in enumerate(CALLERS):
            positions = [
                metric_index + (-0.18 if index == 0 else 0.18)
                for metric_index in range(len(METRICS))
            ]
            values = [
                float(records[caller, variant_type][field])
                for field, _ in METRICS
            ]
            bars = ax.barh(
                positions, values, height=0.30,
                color=color, label=label, zorder=3,
            )
            if variant_type == "SNP":
                handles.append(bars)
            for bar, value in zip(bars, values):
                ax.text(
                    value - 0.018,
                    bar.get_y() + bar.get_height() / 2,
                    f"{value:.2%}",
                    ha="right", va="center",
                    color="white", fontsize=9, fontweight="bold",
                )
        ax.set_yticks(range(len(METRICS)), [label for _, label in METRICS])
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.set_xlabel("Score")
        ax.grid(axis="x", color="#E5E7EB", linewidth=0.8, zorder=0)
        ax.tick_params(axis="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.suptitle(
        "HG002 small-variant benchmarking",
        x=0.09, y=0.97, ha="left", fontsize=18, fontweight="bold",
    )
    fig.text(
        0.09, 0.885,
        "Selected chr20 CMRG gene regions · GIAB v5.0q draft benchmark",
        fontsize=11, color="#4B5563",
    )
    fig.legend(
        handles, [label for _, label, _ in CALLERS],
        loc="upper left", bbox_to_anchor=(0.08, 0.865),
        ncol=2, frameon=False,
    )
    fig.text(
        0.09, 0.115,
        "Evaluable territory: 188,823 bp. Published callsets; supplied filters differ.",
        fontsize=9, color="#4B5563",
    )
    fig.text(
        0.09, 0.065,
        "Descriptive results for one sample and selected regions; "
        "not genome-wide or clinical validation.",
        fontsize=9, color="#4B5563",
    )

    FIGURES.mkdir(parents=True, exist_ok=True)
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    outputs = []
    for extension in ("png", "svg", "pdf"):
        target = FIGURES / f"overall_benchmark_metrics.{extension}"
        fig.savefig(target, dpi=300, facecolor="white")
        outputs.append({
            "path": target.relative_to(ROOT).as_posix(),
            "sha256": sha256(target),
        })
        print("Saved:", target.relative_to(ROOT))
    plt.close(fig)

    metadata = {
        "script": Path(__file__).relative_to(ROOT).as_posix(),
        "script_sha256": sha256(Path(__file__)),
        "inputs": inputs,
        "row_selection": {"Filter": "PASS", "Type": ["SNP", "INDEL"]},
        "metric_columns": [field for field, _ in METRICS],
        "axis_range": [0, 1],
        "outputs": outputs,
        "note": "Additional figure provenance; existing export manifest unchanged.",
    }
    target = PROVENANCE / "overall_benchmark_figure.json"
    target.write_text(json.dumps(metadata, indent=2) + "\n")
    print("Saved:", target.relative_to(ROOT))


if __name__ == "__main__":
    main()
