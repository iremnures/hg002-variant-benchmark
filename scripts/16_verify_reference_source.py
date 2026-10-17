#!/usr/bin/env python3
"""Verify normalized chr20 sequence against a versioned NCBI accession."""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
ACCESSION = "NC_000020.11"
EXPECTED_LENGTH = 64_444_167
SOURCE_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=NC_000020.11&rettype=fasta&retmode=text"
)


def inspect_fasta(handle, label):
    file_hash = hashlib.sha256()
    sequence_hash = hashlib.sha256()
    header = None
    length = 0
    next_progress = 8_000_000
    allowed = set(b"ACGTRYSWKMBDHVN")

    for line in handle:
        file_hash.update(line)
        stripped = line.strip()
        if not stripped:
            continue

        if stripped.startswith(b">"):
            if header is not None:
                raise ValueError(f"{label}: multiple FASTA records")
            header = stripped[1:].decode("ascii")
            continue

        if header is None:
            raise ValueError(f"{label}: missing FASTA header")

        sequence = b"".join(stripped.split()).upper()
        if not set(sequence).issubset(allowed):
            raise ValueError(f"{label}: invalid sequence characters")

        length += len(sequence)
        if length > EXPECTED_LENGTH:
            raise ValueError(f"{label}: sequence exceeds expected length")
        sequence_hash.update(sequence)

        if length >= next_progress:
            print(f"{label}: {length:,} baz okundu", flush=True)
            next_progress += 8_000_000

    if header is None or length == 0:
        raise ValueError(f"{label}: empty FASTA")

    return {
        "header": header,
        "length_bp": length,
        "file_sha256": file_hash.hexdigest(),
        "sequence_sha256": sequence_hash.hexdigest(),
    }


def main():
    manifest_path = ROOT / "config/inputs.json"
    report_path = ROOT / "results/validation/reference_source_verification.json"
    temporary = ROOT / "config/inputs.reference_update.tmp"

    if report_path.exists() or temporary.exists():
        raise FileExistsError("Existing report or temporary update protected")

    original = manifest_path.read_bytes()
    manifest = json.loads(original)
    entries = [
        item for item in manifest["source_inputs"]
        if item["id"] == "reference_chr20"
    ]
    if len(entries) != 1:
        raise ValueError("Expected exactly one chr20 reference input")
    entry = entries[0]
    local_path = ROOT / entry["local_path"]

    before = local_path.stat()
    print("Yerel referans okunuyor...", flush=True)
    with local_path.open("rb") as handle:
        local = inspect_fasta(handle, "Yerel")
    after = local_path.stat()

    if (before.st_size, before.st_mtime_ns) != (
        after.st_size, after.st_mtime_ns
    ):
        raise RuntimeError("Local reference changed during validation")
    if local["file_sha256"] != entry["sha256"]:
        raise ValueError("Local reference differs from input manifest")
    if local["header"].split()[0] != "chr20":
        raise ValueError("Expected local FASTA record name chr20")

    print("\nNCBI chr20 dizisi okunuyor; mevcut dosya değiştirilmiyor...",
          flush=True)
    with urlopen(SOURCE_URL, timeout=30) as response:
        resolved_url = response.geturl()
        remote = inspect_fasta(response, "NCBI")

    matched = (
        remote["header"].split()[0] == ACCESSION
        and local["length_bp"] == remote["length_bp"] == EXPECTED_LENGTH
        and local["sequence_sha256"] == remote["sequence_sha256"]
    )

    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "accession": ACCESSION,
        "source_url": SOURCE_URL,
        "resolved_url": resolved_url,
        "retrieval_note": (
            "Accession-based retrieval URL reconstructed for verification; "
            "the original complete download query was not retained."
        ),
        "normalization": (
            "Ignore FASTA header, whitespace and sequence letter case"
        ),
        "local_path": entry["local_path"],
        "local": local,
        "remote": remote,
        "status": "VERIFIED_SEQUENCE" if matched else "MISMATCH",
        "input_manifest_sha256": hashlib.sha256(original).hexdigest(),
        "script_sha256": hashlib.sha256(
            Path(__file__).resolve().read_bytes()
        ).hexdigest(),
        "scope": "chr20 sequence identity, not whole reference-bundle identity",
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print("\n" + report["status"])
    if not matched:
        raise SystemExit("Mismatch recorded; manifest was not updated.")

    entry.update({
        "source_url": SOURCE_URL,
        "source_status": "remote_accession_sequence_sha256_verified",
        "reference_accession": ACCESSION,
        "sequence_length_bp": local["length_bp"],
        "sequence_sha256": local["sequence_sha256"],
        "normalization": report["normalization"],
        "preparation_note": "Local FASTA record identifier is chr20",
        "provenance_record": str(report_path.relative_to(ROOT)),
    })

    if manifest_path.read_bytes() != original:
        raise RuntimeError("Manifest changed during verification")

    with temporary.open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    temporary.replace(manifest_path)

    verified_states = {
        "remote_size_and_md5_verified",
        "remote_published_md5_verified",
        "pinned_source_sha256_verified",
        "remote_accession_sequence_sha256_verified",
    }
    verified = sum(
        item["source_status"] in verified_states
        for item in manifest["source_inputs"]
    )
    print(f"Saved: {report_path.relative_to(ROOT)}")
    print(f"Updated: config/inputs.json — verified inputs: {verified}/7")


if __name__ == "__main__":
    main()
