from pathlib import Path
import argparse
import json

parser = argparse.ArgumentParser(description="Prepare AD-free comparison copies.")
parser.add_argument("--input-dir", default="data/gene_subset")
parser.add_argument("--output-dir")
args = parser.parse_args()

folder = Path(args.input_dir)
checks = json.loads((folder / "reference_check.json").read_text())
callers = ("giab", "deepvariant", "gatk")

for caller in callers:
    check = checks[caller]
    if check["records_checked"] == 0 or check["mismatches"]:
        raise SystemExit(f"Referans kontrolu gecmemis: {caller}")

outdir = Path(args.output_dir) if args.output_dir else folder / "noAD"
targets = {caller: outdir / f"{caller}.vcf" for caller in callers}
report = outdir / "preparation.json"

if any(p.exists() for p in [*targets.values(), report]):
    raise SystemExit("Hedeflerden biri zaten var; degistirilmedi.")

outdir.mkdir(parents=True, exist_ok=True)
results = {}

for caller in callers:
    source = folder / f"{caller}.padded.vcf"
    count = 0
    changed = 0

    with source.open() as reader, targets[caller].open("x") as writer:
        for line in reader:
            if line.startswith("##FORMAT=<ID=AD,"):
                continue
            if line.startswith("#"):
                writer.write(line)
                continue

            f = line.rstrip("\n").split("\t")
            keys = f[8].split(":")
            if "AD" in keys:
                keep = [i for i, key in enumerate(keys) if key != "AD"]
                f[8] = ":".join(keys[i] for i in keep)
                for column in range(9, len(f)):
                    if f[column] == ".":
                        continue
                    values = f[column].split(":")
                    f[column] = ":".join(
                        values[i] if i < len(values) else "."
                        for i in keep
                    )
                changed += 1

            writer.write("\t".join(f) + "\n")
            count += 1

    if count != checks[caller]["records_checked"]:
        raise SystemExit(f"Kayit sayisi beklenenle uyusmuyor: {caller}")

    results[caller] = {
        "source": str(source),
        "output": str(targets[caller]),
        "records_retained": count,
        "records_with_AD_removed": changed
    }
    print(caller, "| Korunan kayit:", count, "| AD cikarilan:", changed)

with report.open("x") as writer:
    json.dump(results, writer, indent=2)
    writer.write("\n")

print("Karsilastirma kopyalari hazir:", outdir)
