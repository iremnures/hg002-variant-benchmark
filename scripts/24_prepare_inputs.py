#!/usr/bin/env python3
"""Verify existing inputs or retrieve missing inputs from recorded sources."""

import argparse
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
VERIFIED_STATES = {
    "remote_size_and_md5_verified",
    "remote_published_md5_verified",
    "pinned_source_sha256_verified",
    "remote_accession_sequence_sha256_verified",
}


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def project_path(root, relative):
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Expected project-relative path: {relative}")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"Path escapes project: {relative}")
    return path


def write_reference(source, destination, recipe):
    width = recipe["line_bases"]
    newline = recipe["newline"].encode("ascii")
    extra = recipe.get("trailing_blank_lines", 0)
    if width <= 0 or newline not in (b"\n", b"\r\n") or not 0 <= extra <= 10:
        raise ValueError("Invalid reference formatting recipe")

    header = source.readline().strip()
    if not header.startswith(b">"):
        raise ValueError("Remote reference is not FASTA")
    accession = header[1:].split()[0].decode("ascii")
    if accession != recipe["reference_accession"]:
        raise ValueError(f"Unexpected reference accession: {accession}")

    destination.write(b">" + recipe["record_name"].encode("ascii") + newline)
    digest = hashlib.sha256()
    pending = b""
    total = 0
    next_progress = 8_000_000
    allowed = set(b"ACGTRYSWKMBDHVN")

    for line in source:
        if line.startswith(b">"):
            raise ValueError("Multiple reference records returned")
        sequence = b"".join(line.split()).upper()
        if not set(sequence).issubset(allowed):
            raise ValueError("Invalid reference sequence characters")
        total += len(sequence)
        if total > recipe["sequence_length_bp"]:
            raise ValueError("Reference exceeds expected length")
        digest.update(sequence)
        pending += sequence

        while len(pending) >= width:
            destination.write(pending[:width] + newline)
            pending = pending[width:]

        if total >= next_progress:
            print(f"Reference: {total:,} bases read", flush=True)
            next_progress += 8_000_000

    if total != recipe["sequence_length_bp"]:
        raise ValueError("Reference length differs")
    if digest.hexdigest() != recipe["expected_sequence_sha256"]:
        raise ValueError("Reference sequence identity differs")

    if pending:
        destination.write(pending)
        if recipe["final_newline"]:
            destination.write(newline)
    elif not recipe["final_newline"]:
        raise ValueError("Unsupported no-final-newline recipe")

    destination.write(newline * extra)


def install_verified(destination, expected_sha256, writer):
    """Write a temporary file; install without overwriting an existing path."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent, prefix=destination.name + ".",
            suffix=".part", delete=False
        ) as handle:
            temporary = Path(handle.name)
            writer(handle)

        if sha256(temporary) != expected_sha256:
            raise ValueError(f"Prepared file identity differs: {destination}")

        # Same-filesystem hard link: fails rather than overwriting a destination.
        os.link(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def retrieve(url, handle, recipe=None):
    with urlopen(url, timeout=30) as response:
        if recipe:
            write_reference(response, handle, recipe)
            return

        total = 0
        next_progress = 32 * 1024 * 1024
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)
            total += len(chunk)
            if total >= next_progress:
                print(f"Downloaded: {total / 1024 / 1024:.0f} MiB", flush=True)
                next_progress += 32 * 1024 * 1024


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--fetch-missing", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()

    manifest_path = root / "config/inputs.json"
    recipe_path = root / "config/reference_recipe.json"
    manifest = json.loads(manifest_path.read_text())
    recipe = json.loads(recipe_path.read_text())
    manifest_hash = sha256(manifest_path)
    recipe_hash = sha256(recipe_path)

    entries = manifest["source_inputs"]
    if len({item["id"] for item in entries}) != len(entries):
        raise ValueError("Duplicate input IDs")

    reference = next(item for item in entries if item["id"] == "reference_chr20")
    index_entry = next(
        item for item in manifest["derived_inputs"]
        if item["id"] == "reference_index"
    )
    if (
        recipe["expected_fasta_sha256"] != reference["sha256"]
        or recipe["local_fasta_path"] != reference["local_path"]
        or recipe["expected_index_sha256"] != index_entry["sha256"]
        or recipe["local_index_path"] != index_entry["local_path"]
    ):
        raise ValueError("Reference recipe conflicts with input manifest")

    report = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "RUNNING",
        "manifest_sha256": manifest_hash,
        "reference_recipe_sha256": recipe_hash,
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": [],
    }

    try:
        for item in entries:
            if item["source_status"] not in VERIFIED_STATES:
                raise ValueError(f"Unverified source: {item['id']}")
            path = project_path(root, item["local_path"])

            if path.exists():
                if not path.is_file() or sha256(path) != item["sha256"]:
                    raise ValueError(f"Existing input differs: {item['id']}")
                status = "VERIFIED_EXISTING"
            else:
                if not args.fetch_missing:
                    raise FileNotFoundError(
                        f"Missing {item['id']}; use --fetch-missing to retrieve"
                    )
                print(f"Retrieving: {item['id']}", flush=True)
                reference_recipe = recipe if item["id"] == "reference_chr20" else None
                install_verified(
                    path, item["sha256"],
                    lambda handle, item=item, reference_recipe=reference_recipe:
                        retrieve(item["source_url"], handle, reference_recipe),
                )
                status = "PREPARED_AND_VERIFIED"

            report["inputs"].append({
                "id": item["id"],
                "path": item["local_path"],
                "sha256": item["sha256"],
                "status": status,
            })
            print(f"{item['id']}: {status}", flush=True)

        index = project_path(root, index_entry["local_path"])
        if index.exists():
            if not index.is_file() or sha256(index) != index_entry["sha256"]:
                raise ValueError("Existing FASTA index differs")
            index_status = "VERIFIED_EXISTING"
        else:
            if not args.fetch_missing:
                raise FileNotFoundError("FASTA index missing; use --fetch-missing")
            newline_size = len(recipe["newline"].encode("ascii"))
            offset = len((">" + recipe["record_name"]).encode("ascii")) + newline_size
            width = recipe["line_bases"]
            text = (
                f"{recipe['record_name']}\t{recipe['sequence_length_bp']}\t"
                f"{offset}\t{width}\t{width + newline_size}\n"
            ).encode()
            install_verified(
                index, index_entry["sha256"], lambda handle: handle.write(text)
            )
            index_status = "PREPARED_AND_VERIFIED"

        report["reference_index_status"] = index_status
        print(f"reference_index: {index_status}")

        if sha256(manifest_path) != manifest_hash or sha256(recipe_path) != recipe_hash:
            raise RuntimeError("Configuration changed during preparation")
        report["status"] = "PASS"

    except Exception as error:
        report["status"] = "FAIL"
        report["error"] = str(error)

    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    report_dir = root / "results/setup"
    report_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    output = report_dir / f"input_preparation_{stamp}.json"
    with output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")

    print(f"\nSTATUS: {report['status']}")
    print(f"Saved: {output}")
    if report["status"] != "PASS":
        raise SystemExit(report["error"])


if __name__ == "__main__":
    main()
