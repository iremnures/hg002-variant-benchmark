#!/usr/bin/env python3
"""Create a curated, traceable snapshot of publication assets."""

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def sanitize(value):
    if isinstance(value, str):
        return value.replace(str(ROOT), "${PROJECT_ROOT}")
    if isinstance(value, list):
        return [sanitize(item) for item in value]
    if isinstance(value, dict):
        return {key: sanitize(item) for key, item in value.items()}
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default="results/nextflow_regions")
    parser.add_argument("--output-dir", default="reports")
    args = parser.parse_args()

    run = Path(args.run_dir).resolve()
    output = Path(args.output_dir).resolve()
    if output.exists():
        raise FileExistsError(f"Existing report snapshot protected: {output}")

    validation = json.loads((run / "validation.json").read_text())
    if validation.get("status") != "PASS":
        raise ValueError("Selected run has not passed validation")

    clean = ROOT / "results/nextflow_clean_run/validation.json"
    if json.loads(clean.read_text()).get("status") != "PASS":
        raise ValueError("Original uncached workflow validation is not PASS")

    sources = {}
    for caller in ("deepvariant", "gatk"):
        for kind in ("summary", "extended"):
            sources[f"tables/{caller}.{kind}.csv"] = (
                run / f"benchmarks/{caller}/{caller}.{kind}.csv"
            )
            sources[f"tables/{caller}.per_gene.{kind}.csv"] = (
                run / f"per_gene/{caller}_per_gene/{caller}.{kind}.csv"
            )

    for name in ("variants.csv", "summary.csv"):
        sources[f"tables/MYT1.{name}"] = run / f"context_analysis/MYT1/{name}"

    sources["tables/gene_coverage.csv"] = run / "regions/coverage.csv"

    for extension in ("png", "svg", "pdf"):
        name = f"MYT1_indel_recall_by_context.{extension}"
        sources[f"figures/{name}"] = run / f"figures/{name}"

    json_sources = {
        "validation/region_workflow.json": run / "validation.json",
        "validation/uncached_workflow.json": clean,
        "validation/ad_content_preservation.json": (
            ROOT / "results/validation/ad_content_preservation.json"
        ),
        "provenance/callset_sources.json": (
            ROOT / "results/validation/callset_source_verification.json"
        ),
        "provenance/giab_sources.json": (
            ROOT / "results/validation/giab_source_verification.json"
        ),
        "provenance/cmrg_source.json": (
            ROOT / "results/validation/cmrg_source_verification.json"
        ),
        "provenance/reference_source.json": (
            ROOT / "results/validation/reference_source_verification.json"
        ),
        "provenance/tandem_source.json": (
            ROOT / "results/validation/tandem_checksum_discovery.json"
        ),
        "provenance/MYT1_analysis_metadata.json": (
            run / "context_analysis/MYT1/metadata.json"
        ),
    }

    for path in [*sources.values(), *json_sources.values()]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty source asset: {path}")

    # Validate JSON before creating the publication directory.
    parsed = {
        name: json.loads(path.read_text())
        for name, path in json_sources.items()
    }

    warning_lines = []
    for name in ("deepvariant_pilot.log", "gatk_pilot_noAD.log"):
        path = ROOT / "logs" / name
        if not path.is_file():
            raise FileNotFoundError(path)
        matches = [
            (number, line)
            for number, line in enumerate(path.read_text().splitlines(), start=1)
            if "too many AD fields" in line
        ]
        if not matches:
            raise ValueError(f"Expected AD-warning evidence not found: {name}")
        warning_lines.append(f"Source: logs/{name}")
        warning_lines.extend(
            f"Line {number}: {line}" for number, line in matches[:4]
        )
        warning_lines.append("")

    output.mkdir(parents=True)
    assets = []

    for name, source in sources.items():
        destination = output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        assets.append({
            "export": name,
            "source": str(source.relative_to(ROOT)),
            "source_sha256": sha256(source),
            "export_sha256": sha256(destination),
            "transformation": "none",
        })

    for name, source in json_sources.items():
        data = parsed[name]
        transformations = ["Replace absolute project root with ${PROJECT_ROOT}"]
        if name == "provenance/tandem_source.json":
            data = dict(data)
            data.pop("checksum_document", None)
            data.pop("tandem_candidates", None)
            transformations.append(
                "Omit full checksum listing and unrelated candidates; "
                "retain listing hash and exact matches"
            )

        destination = output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x") as handle:
            json.dump(sanitize(data), handle, indent=2)
            handle.write("\n")

        assets.append({
            "export": name,
            "source": str(source.relative_to(ROOT)),
            "source_sha256": sha256(source),
            "export_sha256": sha256(destination),
            "transformation": transformations,
        })

    evidence = output / "provenance/ad_warning_examples.txt"
    evidence.write_text("\n".join(warning_lines))
    assets.append({
        "export": str(evidence.relative_to(output)),
        "source": [
            "logs/deepvariant_pilot.log", "logs/gatk_pilot_noAD.log"
        ],
        "export_sha256": sha256(evidence),
        "transformation": "First four matching warning lines from each source log",
    })

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_run": str(run.relative_to(ROOT)),
        "source_run_task_status_counts": validation.get("task_status_counts", {}),
        "script_sha256": sha256(Path(__file__).resolve()),
        "scope": (
            "Curated snapshot of the selected validated run. Separate earlier "
            "nine-task uncached evidence is also retained."
        ),
        "assets": assets,
    }
    with (output / "manifest.json").open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    print(f"Saved: {output.relative_to(ROOT)}")
    print(f"Exported assets: {len(assets)}")
    print("Original results and logs were not modified.")


if __name__ == "__main__":
    main()
