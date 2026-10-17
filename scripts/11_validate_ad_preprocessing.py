#!/usr/bin/env python3
"""Verify that AD removal preserves all other VCF content."""

import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def comparable_lines(handle):
    for line in handle:
        if not line.startswith("##FORMAT=<ID=AD,"):
            yield line.rstrip("\n")


def sample_values(keys, sample):
    if sample == ".":
        return {key: "." for key in keys}
    values = sample.split(":")
    if len(values) > len(keys):
        raise ValueError("More sample values than FORMAT keys")
    return {
        key: values[i] if i < len(values) else "."
        for i, key in enumerate(keys)
    }


def validate(caller, expected):
    source = ROOT / f"data/gene_subset/{caller}.padded.vcf"
    target = ROOT / f"data/gene_subset/noAD/{caller}.vcf"
    initial_hashes = (sha256(source), sha256(target))
    count, ad_removed = 0, 0
    filters = Counter()

    with source.open() as original, target.open() as prepared:
        pairs = zip_longest(
            comparable_lines(original), comparable_lines(prepared)
        )
        for number, (before, after) in enumerate(pairs, start=1):
            if before is None or after is None:
                raise ValueError(f"{caller}: unequal line counts")
            if before.startswith("#") or after.startswith("#"):
                if before != after:
                    raise ValueError(f"{caller}: header mismatch at line {number}")
                continue

            a, b = before.split("\t"), after.split("\t")
            if len(a) != len(b) or len(a) < 10:
                raise ValueError(f"{caller}: invalid or unequal column counts")
            locus = f"{a[0]}:{a[1]}"
            if a[:8] != b[:8]:
                raise ValueError(f"{caller}: first eight fields changed at {locus}")

            old_keys, new_keys = a[8].split(":"), b[8].split(":")
            if len(set(old_keys)) != len(old_keys):
                raise ValueError(f"{caller}: duplicate FORMAT keys at {locus}")
            if new_keys != [key for key in old_keys if key != "AD"]:
                raise ValueError(f"{caller}: unexpected FORMAT change at {locus}")

            ad_removed += int("AD" in old_keys)
            for old_sample, new_sample in zip(a[9:], b[9:]):
                old = sample_values(old_keys, old_sample)
                new = sample_values(new_keys, new_sample)
                old.pop("AD", None)
                if old != new:
                    raise ValueError(f"{caller}: sample values changed at {locus}")

            count += 1
            filters[a[6]] += 1

    if count != expected["records_retained"]:
        raise ValueError(f"{caller}: count differs from preparation report")
    if ad_removed != expected["records_with_AD_removed"]:
        raise ValueError(f"{caller}: AD count differs from preparation report")
    if initial_hashes != (sha256(source), sha256(target)):
        raise ValueError(f"{caller}: input changed during validation")

    return {
        "status": "PASS",
        "source": str(source.relative_to(ROOT)),
        "output": str(target.relative_to(ROOT)),
        "source_sha256": initial_hashes[0],
        "output_sha256": initial_hashes[1],
        "records_checked": count,
        "records_with_AD_removed": ad_removed,
        "filter_distribution": dict(sorted(filters.items())),
    }


def main():
    report = ROOT / "results/validation/ad_content_preservation.json"
    if report.exists():
        raise FileExistsError(f"Existing report protected: {report}")

    preparation_path = ROOT / "data/gene_subset/noAD/preparation.json"
    preparation_hash = sha256(preparation_path)
    preparation = json.loads(preparation_path.read_text())

    results = {}
    for caller in ("giab", "deepvariant", "gatk"):
        results[caller] = validate(caller, preparation[caller])
        print(f"PASS: {caller} — {results[caller]['records_checked']} records")

    if preparation_hash != sha256(preparation_path):
        raise ValueError("Preparation report changed during validation")

    payload = {
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "preparation_report_sha256": preparation_hash,
        "checks": [
            "Headers equal after excluding FORMAT/AD definition",
            "Record order and first eight VCF columns preserved",
            "Only AD removed from FORMAT keys",
            "Non-AD sample values preserved",
            "Missing trailing FORMAT values treated as missing",
            "Counts agree with preparation report",
            "Input file hashes stable during validation",
        ],
        "scope": "Content preservation; not benchmark outcome equivalence",
        "datasets": results,
    }

    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("x") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    print(f"Saved: {report.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
