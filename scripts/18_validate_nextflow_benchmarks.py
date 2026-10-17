#!/usr/bin/env python3
"""Compare Nextflow benchmarks with the original benchmark outputs."""

import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_summary(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        rows = {}
        for row in reader:
            key = (row["Type"], row["Filter"])
            if key in rows:
                raise ValueError(f"Duplicate summary row: {path}, {key}")
            rows[key] = row
    return columns, rows


def equivalent(a, b, column):
    if a == b:
        return True
    try:
        left, right = Decimal(a), Decimal(b)
    except InvalidOperation:
        return False

    if left.is_nan() or right.is_nan():
        return left.is_nan() and right.is_nan()
    if not left.is_finite() or not right.is_finite():
        return left == right

    if column.startswith("METRIC."):
        return abs(left - right) <= Decimal("1e-9")
    return left == right


def vcf_records(path):
    records = Counter()
    samples = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#") or not line.strip():
                continue
            if samples is None:
                raise ValueError(f"Missing VCF header: {path}")

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 + len(samples):
                raise ValueError(f"Unexpected VCF column count: {path}")

            keys = fields[8].split(":")
            sample_values = {}
            for name, sample in zip(samples, fields[9:]):
                values = sample.split(":")
                if len(values) > len(keys):
                    raise ValueError(f"Invalid sample FORMAT values: {path}")
                sample_values[name] = {
                    key: values[i] if i < len(values) else "."
                    for i, key in enumerate(keys)
                }

            # Header timestamps and command paths are intentionally excluded.
            # INFO and FORMAT field ordering does not affect comparison.
            canonical = {
                "fields": fields[:7],
                "info": sorted(fields[7].split(";")),
                "samples": sample_values,
            }
            records[json.dumps(canonical, sort_keys=True)] += 1
    return records


def main():
    output = ROOT / "results/nextflow_benchmark/validation.json"
    if output.exists():
        raise FileExistsError(f"Existing validation protected: {output}")

    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "status": "PASS",
        "scope": "Summary metrics and canonical VCF record equivalence",
        "callers": {},
    }

    for caller in ("deepvariant", "gatk"):
        baseline = ROOT / f"results/gene_workflow/{caller}"
        produced = ROOT / f"results/nextflow_benchmark/benchmarks/{caller}/{caller}"

        old_csv = Path(str(baseline) + ".summary.csv")
        new_csv = Path(str(produced) + ".summary.csv")
        old_vcf = Path(str(baseline) + ".vcf.gz")
        new_vcf = Path(str(produced) + ".vcf.gz")

        old_columns, old_rows = read_summary(old_csv)
        new_columns, new_rows = read_summary(new_csv)
        differences = []

        if set(old_columns) != set(new_columns):
            differences.append("Summary columns differ")
        if set(old_rows) != set(new_rows):
            differences.append("Summary row identities differ")

        for key in sorted(set(old_rows) & set(new_rows)):
            for column in sorted(set(old_columns) & set(new_columns)):
                a, b = old_rows[key][column], new_rows[key][column]
                if not equivalent(a, b, column):
                    differences.append({
                        "row": key,
                        "column": column,
                        "baseline": a,
                        "nextflow": b,
                    })

        old_records = vcf_records(old_vcf)
        new_records = vcf_records(new_vcf)
        records_equal = old_records == new_records
        passed = not differences and records_equal

        report["callers"][caller] = {
            "status": "PASS" if passed else "FAIL",
            "summary_differences": differences,
            "baseline_vcf_records": sum(old_records.values()),
            "nextflow_vcf_records": sum(new_records.values()),
            "canonical_vcf_records_equal": records_equal,
            "baseline_only_record_count": sum((old_records - new_records).values()),
            "nextflow_only_record_count": sum((new_records - old_records).values()),
            "file_sha256": {
                str(p.relative_to(ROOT)): sha256(p)
                for p in (old_csv, new_csv, old_vcf, new_vcf)
            },
        }
        if not passed:
            report["status"] = "FAIL"

        print(f"{caller}: {report['callers'][caller]['status']}")
        print(f"  Özet metrikler aynı: {not differences}")
        print(f"  VCF kayıt içerikleri aynı: {records_equal}")

    with output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print(f"\nSONUÇ: {report['status']}")
    print(f"Saved: {output.relative_to(ROOT)}")
    if report["status"] != "PASS":
        raise SystemExit("Farklar kaydedildi; sonraki aşamadan önce incelenmeli.")


if __name__ == "__main__":
    main()
