from pathlib import Path
import argparse
import json

parser = argparse.ArgumentParser(description="Validate VCF REF alleles.")
parser.add_argument("--reference", default="data/reference/chr20.GRCh38.fa")
parser.add_argument("--input-dir", default="data/gene_subset")
parser.add_argument("--report")
args = parser.parse_args()

reference = Path(args.reference)
index = Path(str(reference) + ".fai")
input_dir = Path(args.input_dir)

entries = {}
for line in index.read_text().splitlines():
    f = line.split("\t")
    entries[f[0]] = tuple(map(int, f[1:5]))

if "chr20" not in entries:
    raise SystemExit("Referans indeksinde chr20 yok.")

length, offset, line_bases, line_bytes = entries["chr20"]
report = (
    Path(args.report) if args.report
    else input_dir / "reference_check.json"
)
if report.exists():
    raise SystemExit("Kontrol raporu zaten var; degistirilmedi.")

results = {}
with reference.open("rb") as fasta:
    for caller in ("giab", "deepvariant", "gatk"):
        source = input_dir / f"{caller}.padded.vcf"
        checked = 0
        mismatches = []

        with source.open() as reader:
            for line in reader:
                if line.startswith("#"):
                    continue

                f = line.rstrip().split("\t")
                start = int(f[1]) - 1
                ref = f[3].upper()
                end = start + len(ref)

                if f[0] != "chr20" or start < 0 or end > length:
                    raise SystemExit(f"Gecersiz konum: {source}, {f[:2]}")

                first = offset + (start // line_bases) * line_bytes + start % line_bases
                last_pos = end - 1
                last = offset + (last_pos // line_bases) * line_bytes + last_pos % line_bases + 1

                fasta.seek(first)
                actual = fasta.read(last - first).replace(b"\n", b"").replace(b"\r", b"").decode().upper()

                checked += 1
                if actual != ref:
                    mismatches.append({
                        "position": int(f[1]),
                        "vcf_ref": ref,
                        "fasta_ref": actual
                    })

        results[caller] = {
            "records_checked": checked,
            "mismatches": mismatches
        }
        print(caller, "| Kontrol edilen:", checked,
              "| REF uyusmazligi:", len(mismatches))

with report.open("x") as writer:
    json.dump(results, writer, indent=2)
    writer.write("\n")

if any(r["records_checked"] == 0 or r["mismatches"] for r in results.values()):
    raise SystemExit("Kontrol gecmedi; analize devam etmeyin. Rapora bakin.")

print("Tum secili kayitlarin REF alanlari referans DNA ile uyumlu.")
print("Kaydedildi:", report)
