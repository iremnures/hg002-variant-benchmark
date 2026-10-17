#!/usr/bin/env python3
"""Verify CMRG coordinates against a commit-pinned NIST source."""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
COMMIT = "bc2e305556b9e377535f0e8d7c947159c1da2a30"
SOURCE_URL = (
    "https://raw.githubusercontent.com/usnistgov/giab-cmrg-benchmarkset/"
    + COMMIT
    + "/data/cmrg_coords/HG002_GRCh38_CMRG_coordinates.bed"
)


def main():
    manifest_path = ROOT / "config/inputs.json"
    report_path = ROOT / "results/validation/cmrg_source_verification.json"
    temporary = ROOT / "config/inputs.cmrg_update.tmp"

    if report_path.exists() or temporary.exists():
        raise FileExistsError("Existing report or temporary update protected")

    original = manifest_path.read_bytes()
    manifest = json.loads(original)
    matches = [
        item for item in manifest["source_inputs"]
        if item["id"] == "cmrg_coordinates"
    ]
    if len(matches) != 1:
        raise ValueError("Expected exactly one CMRG input")
    entry = matches[0]

    local = (ROOT / entry["local_path"]).read_bytes()
    local_sha256 = hashlib.sha256(local).hexdigest()
    if local_sha256 != entry["sha256"]:
        raise ValueError("Local file differs from recorded input manifest")

    print("Commit'e sabitlenmiş NIST koordinatları okunuyor...", flush=True)
    with urlopen(SOURCE_URL, timeout=20) as response:
        remote = response.read(2_000_001)
        resolved_url = response.geturl()

    if len(remote) > 2_000_000:
        raise ValueError("Unexpected source file size")

    remote_sha256 = hashlib.sha256(remote).hexdigest()
    matched = local == remote
    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": SOURCE_URL,
        "resolved_url": resolved_url,
        "source_commit": COMMIT,
        "local_path": entry["local_path"],
        "local_size_bytes": len(local),
        "remote_size_bytes": len(remote),
        "local_sha256": local_sha256,
        "remote_sha256": remote_sha256,
        "status": "VERIFIED_EXACT_FILE" if matched else "MISMATCH",
        "input_manifest_sha256": hashlib.sha256(original).hexdigest(),
        "script_sha256": hashlib.sha256(
            Path(__file__).resolve().read_bytes()
        ).hexdigest(),
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print(report["status"])
    if not matched:
        raise SystemExit("Mismatch recorded; manifest was not updated.")

    entry.update({
        "source_url": SOURCE_URL,
        "source_status": "pinned_source_sha256_verified",
        "source_commit": COMMIT,
        "source_sha256": remote_sha256,
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
    }
    verified = sum(
        item["source_status"] in verified_states
        for item in manifest["source_inputs"]
    )

    print(f"Saved: {report_path.relative_to(ROOT)}")
    print(f"Updated: config/inputs.json — verified inputs: {verified}/7")


if __name__ == "__main__":
    main()
