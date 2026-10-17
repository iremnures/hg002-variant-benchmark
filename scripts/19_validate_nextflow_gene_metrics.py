#!/usr/bin/env python3
"""Compare original and Nextflow per-gene quantification tables."""

import csv
import hashlib
import json
import runpy
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/18_validate_nextflow_benchmarks.py"
helpers = runpy.run_path(str(HELPER), run_name="validation_helpers")
equivalent = helpers["equivalent"]

DIMENSIONS = (
    "Type", "Subtype", "Subset", "Filter",
    "Genotype", "QQ.Field", "QQ"
)


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_table(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        if not columns:
            raise ValueError(f"Missing CSV header: {path}")
        dimensions = [key for key in DIMENSIONS if key in columns]
        if "Type" not in dimensions or "Filter" not in dimensions:
            raise ValueError(f"Missing table identity fields: {path}")

        rows = {}
        for row in reader:
            key = tuple(row[field] for field in dimensions)
            if key in rows:
                raise ValueError(f"Duplicate row identity in {path}: {key}")
            rows[key] = row
    return columns, dimensions, rows


def main():
    output = ROOT / "results/nextflow_gene_workflow/validation.json"
    if output.exists():
        raise FileExistsError(f"Existing validation protected: {output}")

    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "scope": "All per-gene summary and extended CSV rows",
        "script_sha256": sha256(Path(__file__).resolve()),
        "helper_sha256": sha256(HELPER),
        "tables": [],
    }

    for caller in ("deepvariant", "gatk"):
        for kind in ("summary", "extended"):
            old_path = ROOT / (
                f"results/gene_workflow/per_gene/{caller}.{kind}.csv"
            )
            new_path = ROOT / (
                f"results/nextflow_gene_workflow/per_gene/"
                f"{caller}_per_gene/{caller}.{kind}.csv"
            )

            old_columns, old_dims, old_rows = read_table(old_path)
            new_columns, new_dims, new_rows = read_table(new_path)
            differences = []

            if set(old_columns) != set(new_columns):
                differences.append({"problem": "Different columns"})
            if old_dims != new_dims:
                differences.append({"problem": "Different identity dimensions"})
            if set(old_rows) != set(new_rows):
                differences.append({
                    "problem": "Different row identities",
                    "baseline_only": len(set(old_rows) - set(new_rows)),
                    "nextflow_only": len(set(new_rows) - set(old_rows)),
                })

            for key in sorted(set(old_rows) & set(new_rows)):
                for column in sorted(set(old_columns) & set(new_columns)):
                    a = old_rows[key][column]
                    b = new_rows[key][column]
                    if not equivalent(a, b, column):
                        differences.append({
                            "row": key,
                            "column": column,
                            "baseline": a,
                            "nextflow": b,
                        })

            status = "PASS" if not differences else "FAIL"
            if status != "PASS":
                report["status"] = "FAIL"

            report["tables"].append({
                "caller": caller,
                "table": kind,
                "status": status,
                "baseline": str(old_path.relative_to(ROOT)),
                "nextflow_output": str(new_path.relative_to(ROOT)),
                "baseline_sha256": sha256(old_path),
                "nextflow_sha256": sha256(new_path),
                "baseline_rows": len(old_rows),
                "nextflow_rows": len(new_rows),
                "differences": differences,
            })

            print(
                f"{caller} / {kind}: {status} "
                f"({len(old_rows)} / {len(new_rows)} satır)"
            )
            for difference in differences[:3]:
                print("  Fark:", difference)

    with output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print(f"\nSONUÇ: {report['status']}")
    print(f"Saved: {output.relative_to(ROOT)}")
    if report["status"] != "PASS":
        raise SystemExit("Tablo farkları incelenmeden ilerlenmemeli.")


if __name__ == "__main__":
    main()
