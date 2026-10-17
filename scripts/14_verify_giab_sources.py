#!/usr/bin/env python3
"""Verify GIAB inputs against published MD5 checksums."""

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit, unquote
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]

TARGETS = {
    "giab_vcf": "HG002_GRCh38_v5.0q_smvar.vcf.gz",
    "giab_bed": "HG002_GRCh38_v5.0q_smvar.benchmark.bed",
}


def fingerprints(path):
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    before = path.stat()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(chunk)
            sha256.update(chunk)

    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (
        after.st_size, after.st_mtime_ns
    ):
        raise RuntimeError(f"File changed during hashing: {path}")

    return {
        "md5": md5.hexdigest(),
        "sha256": sha256.hexdigest(),
        "size_bytes": after.st_size,
    }


def basename(url):
    return Path(unquote(urlsplit(url).path)).name


def main():
    report_path = ROOT / "results/validation/giab_source_verification.json"
    manifest_path = ROOT / "config/inputs.json"
    temporary = manifest_path.with_name("inputs.giab_update.tmp")

    if report_path.exists() or temporary.exists():
        raise FileExistsError("Existing report or temporary update protected")

    discovery = json.loads(
        (ROOT / "results/validation/source_directory_discovery.json").read_text()
    )
    links = {
        link
        for directory in discovery["directories"]
        if directory.get("label") == "giab_v5.0q"
        and directory.get("status") == "LISTED"
        for link in directory.get("links", [])
    }

    def find_url(name):
        matches = [url for url in links if basename(url) == name]
        if len(matches) != 1:
            raise ValueError(f"Expected one discovered source URL: {name}")
        return matches[0]

    checksum_url = find_url("checksum.md5")
    original_manifest = manifest_path.read_bytes()
    manifest = json.loads(original_manifest)
    entries = {item["id"]: item for item in manifest["source_inputs"]}

    print("Resmî checksum dosyası okunuyor...", flush=True)
    with urlopen(checksum_url, timeout=20) as response:
        raw = response.read(1_000_001)
        resolved_url = response.geturl()

    if len(raw) > 1_000_000:
        raise ValueError("Unexpected checksum file size")

    published = {}
    for line in raw.decode("utf-8").splitlines():
        match = re.match(r"^\s*([0-9a-fA-F]{32})\s+\*?(.+?)\s*$", line)
        if not match:
            continue
        digest, name = match.groups()
        name = Path(name).name
        if name in published and published[name] != digest.lower():
            raise ValueError(f"Conflicting published checksums: {name}")
        published[name] = digest.lower()

    results = {}
    for identifier, remote_name in TARGETS.items():
        source_url = find_url(remote_name)
        if remote_name not in published:
            raise ValueError(f"No published checksum found: {remote_name}")

        item = entries[identifier]
        path = ROOT / item["local_path"]
        print(f"{identifier}: yerel dosya kontrol ediliyor...", flush=True)
        identity = fingerprints(path)

        if identity["sha256"] != item["sha256"]:
            raise ValueError(f"File changed since input inventory: {identifier}")
        if identity["size_bytes"] != item["size_bytes"]:
            raise ValueError(f"Size changed since input inventory: {identifier}")

        matched = identity["md5"] == published[remote_name]
        results[identifier] = {
            "local_path": item["local_path"],
            "source_url": source_url,
            "published_md5": published[remote_name],
            **identity,
            "status": "VERIFIED_EXACT_FILE" if matched else "MISMATCH",
        }
        print(f"{identifier}: {results[identifier]['status']}", flush=True)

    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "checksum_url": checksum_url,
        "resolved_checksum_url": resolved_url,
        "checksum_document_sha256": hashlib.sha256(raw).hexdigest(),
        "checksum_document": raw.decode("utf-8"),
        "input_manifest_sha256": hashlib.sha256(original_manifest).hexdigest(),
        "script_sha256": fingerprints(Path(__file__).resolve())["sha256"],
        "datasets": results,
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    if not all(r["status"] == "VERIFIED_EXACT_FILE" for r in results.values()):
        raise SystemExit("Mismatch recorded; input manifest was not updated.")

    for identifier, result in results.items():
        entries[identifier].update({
            "source_url": result["source_url"],
            "source_status": "remote_published_md5_verified",
            "published_md5": result["published_md5"],
            "provenance_record": str(report_path.relative_to(ROOT)),
        })

    if manifest_path.read_bytes() != original_manifest:
        raise RuntimeError("Manifest changed during verification; update stopped")

    with temporary.open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    temporary.replace(manifest_path)

    print(f"Saved: {report_path.relative_to(ROOT)}")
    print("Updated: config/inputs.json — GIAB source identities verified")


if __name__ == "__main__":
    main()
