#!/usr/bin/env python3
"""Check local callsets against candidate public GCS object metadata."""

import base64
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
BUCKET = "brain-genomics-public"
PREFIX = "research/sequencing/grch38/vcf/novaseq/wgs_pcr_free/30x/"

CANDIDATES = {
    "deepvariant": (
        "data/callsets/HG002.deepvariant1.0.novaseq30x.grch38.vcf.gz",
        "HG002.novaseq.pcr-free.30x.deepvariant-v1.0.grch38.vcf.gz",
    ),
    "gatk": (
        "data/callsets/HG002.gatk4.novaseq30x.grch38.vcf.gz",
        "HG002.novaseq.pcr-free.30x.gatk4.grch38.vcf.gz",
    ),
}


def file_hashes(path):
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(chunk)
            sha256.update(chunk)
    return {
        "md5_base64": base64.b64encode(md5.digest()).decode(),
        "sha256": sha256.hexdigest(),
    }


def main():
    output = ROOT / "results/validation/callset_source_verification.json"
    if output.exists():
        raise FileExistsError(f"Existing report protected: {output}")

    for local_name, _ in CANDIDATES.values():
        if not (ROOT / local_name).is_file():
            raise FileNotFoundError(ROOT / local_name)

    results = {}
    for caller, (local_name, remote_name) in CANDIDATES.items():
        path = ROOT / local_name
        object_name = PREFIX + remote_name
        metadata_url = (
            f"https://storage.googleapis.com/storage/v1/b/{BUCKET}/o/"
            + quote(object_name, safe="")
        )
        download_url = (
            f"https://storage.googleapis.com/{BUCKET}/"
            + quote(object_name, safe="/")
        )

        print(f"\n{caller}: yerel dosya kimlikleri hesaplanıyor...", flush=True)
        before = path.stat()
        local = file_hashes(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (
            after.st_size, after.st_mtime_ns
        ):
            raise RuntimeError(f"File changed during hashing: {path}")

        entry = {
            "local_path": local_name,
            "local_size_bytes": after.st_size,
            "local_hashes": local,
            "candidate_metadata_url": metadata_url,
            "candidate_download_url": download_url,
        }

        try:
            with urlopen(metadata_url, timeout=30) as response:
                remote = json.load(response)

            entry["remote_metadata"] = {
                key: remote.get(key)
                for key in (
                    "bucket", "name", "size", "md5Hash", "crc32c",
                    "generation", "updated", "timeCreated", "etag"
                )
            }

            size_match = int(remote["size"]) == after.st_size
            md5_match = (
                local["md5_base64"] == remote["md5Hash"]
                if remote.get("md5Hash") else None
            )
            identity_match = (
                remote.get("bucket") == BUCKET
                and remote.get("name") == object_name
            )

            entry["size_match"] = size_match
            entry["md5_match"] = md5_match
            entry["object_identity_match"] = identity_match

            if identity_match and size_match and md5_match is True:
                entry["status"] = "VERIFIED_EXACT_FILE"
            elif md5_match is None:
                entry["status"] = "UNVERIFIED_NO_REMOTE_MD5"
            else:
                entry["status"] = "MISMATCH"

        except (HTTPError, URLError, TimeoutError) as error:
            entry["status"] = "METADATA_LOOKUP_FAILED"
            entry["error"] = str(error)

        results[caller] = entry
        print(f"{caller}: {entry['status']}", flush=True)

    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": file_hashes(Path(__file__).resolve())["sha256"],
        "method": "Local compressed-file size and MD5 versus GCS object metadata",
        "scope": (
            "Exact file identity only. Does not establish model checkpoint, "
            "training independence or identical caller input BAMs."
        ),
        "callsets": results,
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print(f"\nSaved: {output.relative_to(ROOT)}")

    if not all(
        entry["status"] == "VERIFIED_EXACT_FILE"
        for entry in results.values()
    ):
        raise SystemExit(
            "Source verification incomplete; inspect the recorded statuses."
        )


if __name__ == "__main__":
    main()
