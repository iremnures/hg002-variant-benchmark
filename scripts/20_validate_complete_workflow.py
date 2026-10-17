#!/usr/bin/env python3
"""Validate a complete workflow run against established project outputs."""

import argparse
import csv
import hashlib
import json
import runpy
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.image as mpimg
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
H18 = ROOT / "scripts/18_validate_nextflow_benchmarks.py"
H19 = ROOT / "scripts/19_validate_nextflow_gene_metrics.py"
helpers18 = runpy.run_path(str(H18), run_name="benchmark_helpers")
helpers19 = runpy.run_path(str(H19), run_name="table_helpers")


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def compare_tables(old_path, new_path):
    read = helpers19["read_table"]
    equivalent = helpers19["equivalent"]
    old_columns, old_dims, old_rows = read(old_path)
    new_columns, new_dims, new_rows = read(new_path)

    differences = []
    if set(old_columns) != set(new_columns):
        differences.append("Column sets differ")
    if old_dims != new_dims:
        differences.append("Identity dimensions differ")
    if set(old_rows) != set(new_rows):
        differences.append("Row identities differ")

    for key in sorted(set(old_rows) & set(new_rows)):
        for column in set(old_columns) & set(new_columns):
            if not equivalent(
                old_rows[key][column], new_rows[key][column], column
            ):
                differences.append(f"{key}: {column}")

    return not differences, {
        "baseline_rows": len(old_rows),
        "new_rows": len(new_rows),
        "differences": differences[:20],
        "baseline_sha256": sha256(old_path),
        "new_sha256": sha256(new_path),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir", default="results/nextflow_clean_run"
    )
    parser.add_argument("--expected-tasks", type=int, default=9)
    parser.add_argument("--allow-cached", action="store_true")
    args = parser.parse_args()
    output = Path(args.results_dir).resolve()
    report_path = output / "validation.json"

    if report_path.exists():
        raise FileExistsError(f"Existing validation protected: {report_path}")

    checks = []

    def check(label, passed, details=None):
        checks.append({
            "check": label,
            "passed": bool(passed),
            "details": details or {},
        })
        print(f"{'PASS' if passed else 'FAIL'}: {label}", flush=True)

    def same_file(label, baseline, produced):
        old_hash, new_hash = sha256(baseline), sha256(produced)
        check(label, old_hash == new_hash, {
            "baseline_sha256": old_hash,
            "new_sha256": new_hash,
        })

    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "results_dir": str(output),
        "pipeline_sha256_at_validation": sha256(ROOT / "main.nf"),
        "config_sha256_at_validation": sha256(ROOT / "nextflow.config"),
        "validator_sha256": sha256(Path(__file__).resolve()),
        "helper_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in (H18, H19)},
        "checks": checks,
    }

    try:
        for caller in ("giab", "deepvariant", "gatk"):
            same_file(
                f"{caller}: padded VCF",
                ROOT / f"data/gene_subset/{caller}.padded.vcf",
                output / f"prepared/subset/{caller}.padded.vcf",
            )
            same_file(
                f"{caller}: AD-free VCF",
                ROOT / f"data/gene_subset/noAD/{caller}.vcf",
                output / f"comparison_inputs/noAD/{caller}.vcf",
            )

        old_reference = json.loads(
            (ROOT / "data/gene_subset/reference_check.json").read_text()
        )
        new_reference = json.loads(
            (output / "reference_validation/reference_check.json").read_text()
        )
        check("Reference-check results", old_reference == new_reference)

        for caller in ("deepvariant", "gatk"):
            for kind in ("summary", "extended"):
                passed, details = compare_tables(
                    ROOT / f"results/gene_workflow/{caller}.{kind}.csv",
                    output / f"benchmarks/{caller}/{caller}.{kind}.csv",
                )
                check(f"{caller}: benchmark {kind}", passed, details)

                passed, details = compare_tables(
                    ROOT / f"results/gene_workflow/per_gene/{caller}.{kind}.csv",
                    output / f"per_gene/{caller}_per_gene/{caller}.{kind}.csv",
                )
                check(f"{caller}: per-gene {kind}", passed, details)

            old_records = helpers18["vcf_records"](
                ROOT / f"results/gene_workflow/{caller}.vcf.gz"
            )
            new_records = helpers18["vcf_records"](
                output / f"benchmarks/{caller}/{caller}.vcf.gz"
            )
            check(
                f"{caller}: benchmark VCF records",
                old_records == new_records,
                {
                    "baseline_records": sum(old_records.values()),
                    "new_records": sum(new_records.values()),
                },
            )

        for name in ("variants.csv", "summary.csv"):
            same_file(
                f"MYT1: {name}",
                ROOT / f"results/context_analysis/MYT1/{name}",
                output / f"context_analysis/MYT1/{name}",
            )

        for extension in ("png", "svg", "pdf"):
            path = output / f"figures/MYT1_indel_recall_by_context.{extension}"
            check(
                f"Figure exists: {extension}",
                path.is_file() and path.stat().st_size > 0,
            )

        old_image = mpimg.imread(
            ROOT / "results/figures/MYT1_indel_recall_by_context.png"
        )
        new_image = mpimg.imread(
            output / "figures/MYT1_indel_recall_by_context.png"
        )
        check("PNG pixel equality", np.array_equal(old_image, new_image))

        region_dir = output / "regions"
        if region_dir.exists():
            region_helper_path = ROOT / "scripts/21_prepare_gene_regions.py"
            region_helpers = runpy.run_path(
                str(region_helper_path), run_name="region_helpers"
            )
            bed_rows = region_helpers["bed_rows"]
            normalized_bed = region_helpers["normalized_bed"]
            baseline_regions = ROOT / "data/regions"

            def gene_rows(path):
                return sorted(
                    (chrom, start, end, tuple(extra))
                    for chrom, start, end, extra in bed_rows(path)
                )

            old_genes = gene_rows(baseline_regions / "chr20_cmrg_genes.bed")
            new_genes = gene_rows(region_dir / "chr20_cmrg_genes.bed")
            check("Generated gene selection and coordinates", old_genes == new_genes)

            check(
                "Generated evaluable interval union",
                normalized_bed(
                    baseline_regions / "chr20_cmrg_evaluable.bed", "chr20"
                ) == normalized_bed(
                    region_dir / "chr20_cmrg_evaluable.bed", "chr20"
                ),
            )

            expected_beds = set()
            for chrom, start, end, extra in old_genes:
                name = extra[0]
                old_bed = baseline_regions / f"per_gene/{name}.bed"
                new_bed = region_dir / f"per_gene/{name}.bed"

                if old_bed.exists():
                    expected_beds.add(old_bed.name)
                    passed = new_bed.exists() and normalized_bed(
                        old_bed, "chr20"
                    ) == normalized_bed(new_bed, "chr20")
                else:
                    passed = not new_bed.exists()
                check(f"Generated region: {name}", passed)

            actual_beds = {p.name for p in (region_dir / "per_gene").glob("*.bed")}
            check("No missing or unexpected per-gene BEDs", actual_beds == expected_beds)
            check(
                "Generated gene stratification list",
                sorted((baseline_regions / "per_gene/genes.tsv").read_text().splitlines())
                == sorted((region_dir / "per_gene/genes.tsv").read_text().splitlines()),
            )
            report["region_helper_sha256"] = sha256(region_helper_path)

        with (output / "execution/trace.tsv").open(newline="") as handle:
            tasks = list(csv.DictReader(handle, delimiter="\t"))

        def task_ok(row):
            if row.get("status") == "COMPLETED":
                return row.get("exit") == "0"
            if args.allow_cached and row.get("status") == "CACHED":
                return row.get("exit") in {"0", "-"}
            return False

        trace_passed = len(tasks) == args.expected_tasks and all(
            task_ok(row) for row in tasks
        )
        check(
            f"Expected {args.expected_tasks} successful task records",
            trace_passed,
            {"allow_cached": args.allow_cached},
        )
        report["task_status_counts"] = {
            status: sum(row.get("status") == status for row in tasks)
            for status in sorted({row.get("status", "") for row in tasks})
        }
        report["task_trace"] = tasks

        print("\n--- Nextflow görev kayıtları ---")
        print("CACHED görevlerdeki süreler bu çalışmanın yeni ölçümü değildir.")
        for row in tasks:
            print(
                f"{row.get('name')} | süre={row.get('realtime')} | "
                f"peak_rss={row.get('peak_rss')}"
            )

    except Exception as error:
        check("Validation completed without exception", False, {
            "error": str(error)
        })

    report["status"] = "PASS" if checks and all(
        item["passed"] for item in checks
    ) else "FAIL"

    with report_path.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print(f"\nSONUÇ: {report['status']}")
    print(f"Saved: {report_path}")
    if report["status"] != "PASS":
        raise SystemExit("Validation differences require review.")


if __name__ == "__main__":
    main()
