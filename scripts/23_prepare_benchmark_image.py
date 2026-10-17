#!/usr/bin/env python3
"""Check the pinned benchmark image; optionally pull it if absent."""

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PODMAN = ["podman", "--storage-opt", "overlay.ignore_chown_errors=true"]


def normalize_id(value):
    return value.removeprefix("sha256:")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pull-if-missing", action="store_true")
    args = parser.parse_args()

    config = json.loads(
        (ROOT / "config/benchmark_container.json").read_text()
    )
    image = config["image"]
    expected = normalize_id(config["expected_image_id"])

    inspected = subprocess.run(
        PODMAN + ["image", "inspect", image],
        capture_output=True, text=True,
    )
    if inspected.returncode != 0:
        if not args.pull_if_missing:
            raise SystemExit(
                "Pinned image not available locally. "
                "Use --pull-if-missing to retrieve the recorded image."
            )
        subprocess.run(PODMAN + ["pull", image], check=True)
        inspected = subprocess.run(
            PODMAN + ["image", "inspect", image],
            capture_output=True, text=True, check=True,
        )

    actual = json.loads(inspected.stdout)[0]
    if normalize_id(actual["Id"]) != expected:
        raise SystemExit("Image ID differs from recorded benchmark image.")

    for field, key in (("Architecture", "architecture"), ("Os", "os")):
        if config.get(key) and actual.get(field) != config[key]:
            raise SystemExit(f"Image platform differs: {field}")

    record = ROOT / "docs/software/happy_image_inspect.json"
    if record.exists():
        previous = json.loads(record.read_text())[0]
        if normalize_id(previous["Id"]) != expected:
            raise SystemExit("Existing workflow image record conflicts.")
        print("Existing workflow image record preserved.")
    else:
        record.parent.mkdir(parents=True, exist_ok=True)
        minimal = {
            key: actual.get(key)
            for key in ("Id", "RepoTags", "RepoDigests", "Architecture", "Os")
        }
        with record.open("x") as handle:
            json.dump([minimal], handle, indent=2)
            handle.write("\n")
        print("Created workflow image record.")

    print("PASS: Pinned image identity and platform match.")
    print("Image:", image)


if __name__ == "__main__":
    main()
