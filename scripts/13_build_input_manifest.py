#!/usr/bin/env python3
"""Inventory project inputs without inventing missing source provenance."""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INPUTS = [
    ("giab_vcf", "benchmark",
     "data/benchmark/HG002_GRCh38_v5.0q_smvar.vcf.gz"),
    ("giab_bed", "benchmark_regions",
     "data/benchmark/HG002_GRCh38_v5.0q_smvar.benchmark.bed"),
    ("reference_chr20", "reference_sequence",
     "data/reference/chr20.GRCh38.fa"),
    ("cmrg_coordinates", "region_annotation",
     "data/annotations/HG002_GRCh38_CMRG_coordinates.bed"),
    ("tandem_repeats", "context_annotation",
     "data/stratifications/GRCh38_AllTandemRepeats.bed.gz"),
]

DERIVED = [
    ("selected_regions", "data/regions/chr20_cmrg_genes.bed"),
    ("evaluable_regions", "data/regions/chr20_cmrg_evaluable.bed"),
    ("gene_stratification", "data/regions/per_gene/genes.tsv"),
    ("reference_index", "data/reference/chr20.GRCh38.fa.fai"),
]


def fingerprint(path):
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (
        after.st_size, after.st_mtime_ns
    ):
        raise RuntimeError(f"File changed during hashing: {path}")
    return {"size_bytes": after.st_size, "sha256": digest.hexdigest()}


def main():
    output = ROOT / "config/inputs.json"
    if output.exists():
        raise FileExistsError(f"Existing manifest protected: {output}")

    callset_path = ROOT / "config/callsets.json"
    callsets = json.loads(callset_path.read_text())
    required = [ROOT / p for _, _, p in INPUTS]
    required += [ROOT / p for _, p in DERIVED]
    required += [ROOT / c["local_path"] for c in callsets["callsets"]]

    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")

    sources = []
    for caller in callsets["callsets"]:
        identity = fingerprint(ROOT / caller["local_path"])
        if identity["sha256"] != caller["sha256"]:
            raise ValueError(f"Callset changed: {caller['id']}")
        if identity["size_bytes"] != caller["size_bytes"]:
            raise ValueError(f"Callset size changed: {caller['id']}")

        sources.append({
            "id": caller["id"],
            "role": "query_callset",
            "local_path": caller["local_path"],
            **identity,
            "source_url": caller["source_url"],
            "source_status": "remote_size_and_md5_verified",
            "provenance_record": "config/callsets.json",
        })
        print(f"{caller['id']}: verified source identity retained")

    for identifier, role, name in INPUTS:
        sources.append({
            "id": identifier,
            "role": role,
            "local_path": name,
            **fingerprint(ROOT / name),
            "source_url": None,
            "source_status": "local_identity_recorded_source_unverified",
        })
        print(f"{identifier}: local identity recorded; source URL pending")

    derived = []
    for identifier, name in DERIVED:
        derived.append({
            "id": identifier,
            "local_path": name,
            **fingerprint(ROOT / name),
            "status": "derived_local_file",
        })

    for path in sorted((ROOT / "data/regions/per_gene").glob("*.bed")):
        derived.append({
            "id": f"gene_region_{path.stem}",
            "local_path": str(path.relative_to(ROOT)),
            **fingerprint(path),
            "status": "derived_local_file",
        })

    report = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": fingerprint(Path(__file__).resolve())["sha256"],
        "callset_manifest_sha256": fingerprint(callset_path)["sha256"],
        "source_inputs": sources,
        "derived_inputs": derived,
        "notes": [
            "A local SHA-256 identifies bytes; it does not verify origin.",
            "Null source URLs are unresolved, not verified download recipes.",
            "Derived-file identities do not establish their generation recipe.",
            "This inventories current inputs, not a complete clean rebuild.",
        ],
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print(f"\nSaved: {output.relative_to(ROOT)}")
    print(f"Source inputs: {len(sources)}")
    print(f"Derived inputs: {len(derived)}")


if __name__ == "__main__":
    main()
