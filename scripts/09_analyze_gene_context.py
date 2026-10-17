#!/usr/bin/env python3
"""Compare truth-side INDEL recall by length and tandem-repeat overlap."""

import argparse
import csv
import gzip
import hashlib
import json
import platform
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LENGTHS = ("1-5", "6-15", "16+")
CONTEXTS = ("repeat_overlap", "no_repeat_overlap")


def bed_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval in {path}: {line.strip()}")
            yield chrom, start, end


def overlaps(chrom, start, end, intervals):
    return any(
        chrom == c and start < right and end > left
        for c, left, right in intervals
    )


def length_group(length):
    return "1-5" if length <= 5 else "6-15" if length <= 15 else "16+"


def read_truth(path, regions, repeats):
    records = {}
    truth_column = None

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                truth_column = 9 + samples.index("TRUTH")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if truth_column is None:
                raise ValueError(f"Missing VCF sample header: {path}")

            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref = fields[0], int(fields[1]), fields[3]
            start, end = pos - 1, pos - 1 + len(ref)

            if not overlaps(chrom, start, end, regions):
                continue

            truth = dict(zip(
                fields[8].split(":"),
                fields[truth_column].split(":")
            ))
            decision = truth.get("BD")
            if truth.get("BVT") != "INDEL" or decision not in {"TP", "FN"}:
                continue

            alleles = [ref] + fields[4].split(",")
            gt = truth.get("GT", ".")
            decoded, changes = [], []

            for index in gt.replace("|", "/").split("/"):
                if index == ".":
                    decoded.append("<MISSING>")
                    continue
                if not index.isdigit() or int(index) >= len(alleles):
                    raise ValueError(f"Invalid truth GT at {chrom}:{pos}: {gt}")

                allele = alleles[int(index)]
                decoded.append(allele)
                if allele == "*" or allele.startswith("<"):
                    continue

                difference = abs(len(allele) - len(ref))
                if difference:
                    changes.append(difference)

            if not changes:
                raise ValueError(f"No measurable INDEL allele at {chrom}:{pos}")

            genotype = tuple(sorted(decoded))
            key = (chrom, pos, ref, genotype)
            if key in records:
                raise ValueError(f"Duplicate truth identity: {key}")

            context = (
                "repeat_overlap"
                if overlaps(chrom, start, end, repeats)
                else "no_repeat_overlap"
            )
            records[key] = {
                "decision": decision,
                "max_length": max(changes),
                "length_group": length_group(max(changes)),
                "context": context,
                "partial_gt": "<MISSING>" in genotype,
                "original_gt": gt,
            }

    return records


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", default="MYT1")
    parser.add_argument("--deepvariant")
    parser.add_argument("--gatk")
    parser.add_argument("--gene-bed")
    parser.add_argument("--repeat-bed")
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    gene = args.gene.upper()

    if not gene or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in gene):
        raise ValueError("Invalid gene name")

    paths = {
        "gene_bed": Path(args.gene_bed).resolve() if args.gene_bed else ROOT / f"data/regions/per_gene/{gene}.bed",
        "repeats": Path(args.repeat_bed).resolve() if args.repeat_bed else ROOT / "data/stratifications/GRCh38_AllTandemRepeats.bed.gz",
        "deepvariant": Path(args.deepvariant).resolve() if args.deepvariant else ROOT / "results/gene_workflow/deepvariant.vcf.gz",
        "gatk": Path(args.gatk).resolve() if args.gatk else ROOT / "results/gene_workflow/gatk.vcf.gz",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    output = Path(args.output_dir).resolve() if args.output_dir else ROOT / f"results/context_analysis/{gene}"
    if output.exists():
        raise FileExistsError(
            f"Output already exists; nothing overwritten: {output}"
        )

    regions = list(bed_rows(paths["gene_bed"]))
    if not regions:
        raise ValueError("Gene BED contains no intervals")

    print(f"{gene}: reading tandem-repeat intervals...", flush=True)
    repeats = [
        row for row in bed_rows(paths["repeats"])
        if overlaps(*row, regions)
    ]
    dv = read_truth(paths["deepvariant"], regions, repeats)
    gatk = read_truth(paths["gatk"], regions, repeats)

    if set(dv) != set(gatk):
        difference = sorted(set(dv) ^ set(gatk))
        raise ValueError(f"Truth identities differ; examples: {difference[:5]}")
    if not dv:
        raise ValueError("No evaluable truth INDEL records found")

    variants = []
    for key in sorted(dv):
        chrom, pos, ref, genotype = key
        d, g = dv[key], gatk[key]
        for field in ("max_length", "length_group", "context", "partial_gt"):
            if d[field] != g[field]:
                raise ValueError(f"Context disagreement at {key}: {field}")

        variants.append({
            "gene": gene,
            "chrom": chrom,
            "pos_1based": pos,
            "ref": ref,
            "truth_alleles_json": json.dumps(genotype),
            "partial_truth_gt": d["partial_gt"],
            "deepvariant_truth_gt": d["original_gt"],
            "gatk_truth_gt": g["original_gt"],
            "max_indel_length_bp": d["max_length"],
            "length_group": d["length_group"],
            "repeat_context": d["context"],
            "deepvariant_BD": d["decision"],
            "gatk_BD": g["decision"],
        })

    summary = []
    for length in LENGTHS:
        for context in CONTEXTS:
            group = [
                row for row in variants
                if row["length_group"] == length
                and row["repeat_context"] == context
            ]
            for caller in ("deepvariant", "gatk"):
                counts = Counter(row[f"{caller}_BD"] for row in group)
                total = len(group)
                summary.append({
                    "gene": gene,
                    "caller": caller,
                    "length_group": length,
                    "repeat_context": context,
                    "truth_total": total,
                    "truth_tp": counts["TP"],
                    "truth_fn": counts["FN"],
                    "recall": counts["TP"] / total if total else "",
                })

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "gene": gene,
        "python": platform.python_version(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "truth_record_count": len(variants),
        "inputs": {
            name: {
                "path": str(path),
                "sha256": sha256(path),
            }
            for name, path in paths.items()
        },
        "methods": {
            "selection": "TRUTH BVT=INDEL and BD in TP,FN; REF overlaps gene BED",
            "identity": "CHROM,POS,REF,sorted decoded TRUTH genotype alleles",
            "phase": "Ignored only for cross-output truth-record identity",
            "missing_gt": "Preserved as <MISSING>; no genotype imputation",
            "length": "Maximum absolute length change among called literal alleles",
            "symbolic_alleles": "Excluded from length calculation",
            "repeat_overlap": "Any REF-span overlap; insertions use anchor base",
            "unit": "Truth VCF record, not an independent biological event",
            "scope": "Exploratory recall analysis; no FP/precision calculation",
        },
    }

    output.mkdir(parents=True)
    write_csv(output / "variants.csv", variants)
    write_csv(output / "summary.csv", summary)
    with (output / "metadata.json").open("x") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")

    for caller in ("deepvariant", "gatk"):
        counts = Counter(row[f"{caller}_BD"] for row in variants)
        print(f"{caller}: {len(variants)} records, "
              f"{counts['TP']} TP, {counts['FN']} FN")
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
