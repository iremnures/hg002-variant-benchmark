#!/usr/bin/env python3
"""Check that parameterized preparation reproduces existing outputs."""

import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CALLERS = ("giab", "deepvariant", "gatk")


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load(path):
    return json.loads(path.read_text())


def main():
    workspace = ROOT / "results/reproducibility/preparation_cli"
    if workspace.exists():
        raise FileExistsError(f"Existing check protected: {workspace}")

    baseline = ROOT / "data/gene_subset"
    relative_files = [
        Path(f"{caller}.padded.vcf") for caller in CALLERS
    ] + [
        Path(f"noAD/{caller}.vcf") for caller in CALLERS
    ]
    report_files = [
        Path("preparation.json"),
        Path("reference_check.json"),
        Path("noAD/preparation.json"),
    ]

    # Önce mevcut doğrulama dayanaklarının kimliklerini kaydet.
    baseline_hashes = {
        str(name): sha256(baseline / name)
        for name in relative_files + report_files
    }

    commands = [
        [
            sys.executable, "scripts/04_prepare_gene_vcfs.py",
            "--output-dir", str(workspace),
        ],
        [
            sys.executable, "scripts/05_check_gene_reference.py",
            "--input-dir", str(workspace),
        ],
        [
            sys.executable, "scripts/06_prepare_gene_benchmark_inputs.py",
            "--input-dir", str(workspace),
        ],
    ]
    script_paths = [
        ROOT / "scripts/04_prepare_gene_vcfs.py",
        ROOT / "scripts/05_check_gene_reference.py",
        ROOT / "scripts/06_prepare_gene_benchmark_inputs.py",
    ]

    workspace.mkdir(parents=True)
    result = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Preparation equivalence; benchmark not rerun",
        "status": "RUNNING",
        "script_sha256": {
            str(p.relative_to(ROOT)): sha256(p) for p in script_paths
        },
        "baseline_sha256": baseline_hashes,
        "steps": [],
        "file_comparisons": [],
    }

    try:
        for number, command in enumerate(commands, start=1):
            print(f"\n--- Hazırlama adımı {number}/3 ---", flush=True)
            started = time.perf_counter()
            completed = subprocess.run(
                command, cwd=ROOT, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True
            )
            elapsed = time.perf_counter() - started

            log = workspace / f"step_{number}.log"
            log.write_text(completed.stdout)
            print(completed.stdout, end="", flush=True)

            result["steps"].append({
                "command": command,
                "returncode": completed.returncode,
                "elapsed_seconds": elapsed,
                "log": str(log.relative_to(ROOT)),
            })
            completed.check_returncode()

        for name in relative_files:
            expected = baseline_hashes[str(name)]
            actual = sha256(workspace / name)
            matched = actual == expected
            result["file_comparisons"].append({
                "file": str(name),
                "baseline_sha256": expected,
                "new_sha256": actual,
                "identical": matched,
            })
            if not matched:
                raise ValueError(f"VCF differs from baseline: {name}")
            print(f"AYNI: {name}")

        # Raporlardaki çıktı yolları değişebilir; bilimsel içeriği karşılaştır.
        old = load(baseline / "preparation.json")
        new = load(workspace / "preparation.json")
        for key in (
            "chromosome", "padding_bp", "windows_bed_coordinates",
            "sources", "record_counts", "filtering_applied"
        ):
            if old[key] != new[key]:
                raise ValueError(f"Preparation report differs: {key}")

        if load(baseline / "reference_check.json") != load(
            workspace / "reference_check.json"
        ):
            raise ValueError("Reference-check results differ")

        old = load(baseline / "noAD/preparation.json")
        new = load(workspace / "noAD/preparation.json")
        for caller in CALLERS:
            for key in ("records_retained", "records_with_AD_removed"):
                if old[caller][key] != new[caller][key]:
                    raise ValueError(f"AD preparation differs: {caller}, {key}")

        for name, expected in baseline_hashes.items():
            if sha256(baseline / name) != expected:
                raise ValueError(f"Baseline changed during check: {name}")

        for path in script_paths:
            if sha256(path) != result["script_sha256"][str(path.relative_to(ROOT))]:
                raise ValueError(f"Script changed during check: {path}")

        result["status"] = "PASS"
        result["report_comparisons"] = "Equivalent except output-path metadata"
        print("\nPASS: Altı VCF birebir aynı; kontrol raporları uyumlu.")
        print("Mevcut dayanak dosyaları değişmedi.")

    except Exception as error:
        result["status"] = "FAIL"
        result["error"] = str(error)
    finally:
        result["finished_utc"] = datetime.now(timezone.utc).isoformat()
        output = workspace / "validation.json"
        with output.open("x") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        print(f"Saved: {output.relative_to(ROOT)}")

    if result["status"] != "PASS":
        raise SystemExit(f"Kontrol geçmedi: {result.get('error')}")


if __name__ == "__main__":
    main()
