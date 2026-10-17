#!/usr/bin/env python3
"""Build chromosome-specific gene regions and benchmark intersections."""

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def bed_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"Invalid BED row: {path}:{number}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid interval: {path}:{number}")
            yield fields[0], start, end, fields[3:]


def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def bed_text(chromosome, intervals):
    return "".join(f"{chromosome}\t{s}\t{e}\n" for s, e in intervals)


def normalized_bed(path, chromosome):
    rows = list(bed_rows(path))
    if any(c != chromosome for c, _, _, _ in rows):
        raise ValueError(f"Unexpected chromosome in {path}")
    return merge((s, e) for _, s, e, _ in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinates", default=str(
        ROOT / "data/annotations/HG002_GRCh38_CMRG_coordinates.bed"))
    parser.add_argument("--benchmark-bed", default=str(
        ROOT / "data/benchmark/HG002_GRCh38_v5.0q_smvar.benchmark.bed"))
    parser.add_argument("--reference-index", default=str(
        ROOT / "data/reference/chr20.GRCh38.fa.fai"))
    parser.add_argument("--chromosome", default="chr20")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--validate-against")
    args = parser.parse_args()

    coordinates = Path(args.coordinates).resolve()
    benchmark = Path(args.benchmark_bed).resolve()
    reference_index = Path(args.reference_index).resolve()
    output = Path(args.output_dir).resolve()
    chromosome = args.chromosome

    if output.exists():
        raise FileExistsError(f"Existing output protected: {output}")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", chromosome):
        raise ValueError("Invalid chromosome identifier")

    lengths = {}
    for line in reference_index.read_text().splitlines():
        fields = line.split("\t")
        if fields[0] in lengths:
            raise ValueError("Duplicate reference-index entry")
        lengths[fields[0]] = int(fields[1])
    if chromosome not in lengths:
        raise ValueError("Selected chromosome missing from reference index")

    input_hashes = {
        str(path): sha256(path)
        for path in (coordinates, benchmark, reference_index)
    }

    genes = []
    names = set()
    for chrom, start, end, extra in bed_rows(coordinates):
        if chrom != chromosome:
            continue
        if len(extra) != 1:
            raise ValueError("Expected exactly one gene-name column")
        name = extra[0]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            raise ValueError(f"Invalid gene identifier: {name}")
        if name in names:
            raise ValueError(f"Duplicate gene identifier: {name}")
        if end > lengths[chromosome]:
            raise ValueError(f"Gene exceeds chromosome length: {name}")
        names.add(name)
        genes.append((start, end, name))

    genes.sort()
    if not genes:
        raise ValueError("No genes found on selected chromosome")

    confident = merge(
        (start, end)
        for chrom, start, end, _ in bed_rows(benchmark)
        if chrom == chromosome
    )
    if any(end > lengths[chromosome] for _, end in confident):
        raise ValueError("Benchmark interval exceeds chromosome length")

    per_gene = {}
    coverage = []
    for start, end, name in genes:
        intervals = merge(
            (max(start, left), min(end, right))
            for left, right in confident
            if start < right and end > left
        )
        per_gene[name] = intervals
        covered = sum(right - left for left, right in intervals)
        coverage.append({
            "gene": name,
            "chromosome": chromosome,
            "source_span_bp": end - start,
            "evaluable_bp": covered,
            "coverage_fraction": covered / (end - start),
            "status": "evaluable_region" if covered else "no_evaluable_region",
        })

    evaluable = merge(
        interval for intervals in per_gene.values() for interval in intervals
    )

    files = {
        f"{chromosome}_cmrg_genes.bed": "".join(
            f"{chromosome}\t{s}\t{e}\t{name}\n"
            for s, e, name in genes
        ),
        f"{chromosome}_cmrg_evaluable.bed": bed_text(chromosome, evaluable),
        "per_gene/genes.tsv": "".join(
            f"{name}\t{name}.bed\n"
            for _, _, name in genes if per_gene[name]
        ),
    }
    for _, _, name in genes:
        if per_gene[name]:
            files[f"per_gene/{name}.bed"] = bed_text(chromosome, per_gene[name])

    # Compute first; only then create new output files.
    output.mkdir(parents=True)
    for relative, text in files.items():
        path = output / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x") as handle:
            handle.write(text)

    with (output / "coverage.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(coverage[0]))
        writer.writeheader()
        writer.writerows(coverage)

    checks = []
    if args.validate_against:
        baseline = Path(args.validate_against).resolve()
        old_genes = [
            (start, end, extra[0])
            for chrom, start, end, extra in bed_rows(
                baseline / f"{chromosome}_cmrg_genes.bed"
            )
            if chrom == chromosome
        ]
        checks.append({
            "check": "Selected genes and coordinates",
            "passed": sorted(old_genes) == genes,
        })
        checks.append({
            "check": "Evaluable interval union",
            "passed": normalized_bed(
                baseline / f"{chromosome}_cmrg_evaluable.bed", chromosome
            ) == evaluable,
        })

        for _, _, name in genes:
            old_path = baseline / f"per_gene/{name}.bed"
            same = (
                normalized_bed(old_path, chromosome) == per_gene[name]
                if old_path.exists()
                else not per_gene[name]
            )
            checks.append({"check": f"{name}: evaluable territory", "passed": same})

        old_tsv = (baseline / "per_gene/genes.tsv").read_text().splitlines()
        checks.append({
            "check": "Gene stratification names and relative paths",
            "passed": sorted(old_tsv) == sorted(files["per_gene/genes.tsv"].splitlines()),
        })

    for name, digest in input_hashes.items():
        if sha256(Path(name)) != digest:
            raise RuntimeError(f"Input changed during generation: {name}")

    passed = all(item["passed"] for item in checks)
    report = {
        "status": ("PASS" if passed else "FAIL") if checks else "GENERATED",
        "chromosome": chromosome,
        "coordinate_system": "0-based, half-open",
        "selection": "All source-coordinate genes on selected chromosome",
        "intersection_policy": "Interval intersection; merge overlapping or touching intervals",
        "genes_selected": len(genes),
        "genes_with_evaluable_regions": sum(bool(v) for v in per_gene.values()),
        "evaluable_union_bp": sum(end - start for start, end in evaluable),
        "input_sha256": input_hashes,
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {
            name: sha256(output / name) for name in [*files, "coverage.csv"]
        },
        "baseline_checks": checks,
    }
    with (output / "region_preparation.json").open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    for item in checks:
        print(f"{'PASS' if item['passed'] else 'FAIL'}: {item['check']}")
    print(f"\nSelected genes: {len(genes)}")
    print(f"Genes with evaluable regions: {report['genes_with_evaluable_regions']}")
    print(f"Evaluable union: {report['evaluable_union_bp']} bp")
    print(f"STATUS: {report['status']}")
    print(f"Saved: {output}")

    if checks and not passed:
        raise SystemExit("Generated regions differ from baseline.")


if __name__ == "__main__":
    main()
