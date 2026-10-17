from pathlib import Path
import csv
import json

folder = Path("results/pilot_workflow/stratified")
metrics = [
    "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
    "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK",
    "METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score"
]
regions = {"TANDEM_REPEAT", "OUTSIDE_TANDEM_REPEAT"}

def read_contexts(path):
    with path.open(newline="") as reader:
        rows = [
            row for row in csv.DictReader(reader)
            if row["Filter"] == "PASS"
            and row["Subtype"] == "*"
            and row["Type"] in {"SNP", "INDEL"}
            and row["Subset"] in regions
        ]
    expected = {(region, kind) for region in regions for kind in ("SNP", "INDEL")}
    result = {(row["Subset"], row["Type"]): row for row in rows}
    if len(rows) != 4 or set(result) != expected:
        raise SystemExit(f"Beklenmeyen tablo yapisi: {path}")
    return result

checks = []
for caller in ("deepvariant", "gatk"):
    old_path = Path(f"results/pilot/{caller}_stratified.extended.csv")
    new_path = folder / f"{caller}.extended.csv"
    old = read_contexts(old_path)
    new = read_contexts(new_path)

    for region, kind in sorted(old):
        for metric in metrics:
            before = float(old[(region, kind)][metric])
            after = float(new[(region, kind)][metric])
            if before != after:
                raise SystemExit(f"FARK: {caller} {region} {kind} {metric}")

        print(caller, region, kind, ": kontrol edilen degerler ayni.")
        checks.append({
            "caller": caller,
            "region": region,
            "variant_type": kind,
            "original": str(old_path),
            "workflow": str(new_path),
            "metrics_checked": metrics,
            "passed": True
        })

target = folder / "reproduction_check.json"
with target.open("x") as writer:
    json.dump(checks, writer, indent=2)
    writer.write("\n")

print("Kontrol kaydedildi:", target)
