from pathlib import Path
import argparse
import gzip
import json

parser = argparse.ArgumentParser(description="Prepare padded gene-region VCFs.")
parser.add_argument("--regions", default="data/regions/chr20_cmrg_genes.bed")
parser.add_argument("--giab", default="data/benchmark/HG002_GRCh38_v5.0q_smvar.vcf.gz")
parser.add_argument("--deepvariant", default="data/callsets/HG002.deepvariant1.0.novaseq30x.grch38.vcf.gz")
parser.add_argument("--gatk", default="data/callsets/HG002.gatk4.novaseq30x.grch38.vcf.gz")
parser.add_argument("--output-dir", default="data/gene_subset")
parser.add_argument("--padding", type=int, default=10_000)
args = parser.parse_args()
if args.padding < 0:
    parser.error("--padding must be nonnegative")

padding = args.padding
chromosome_length = 64_444_167
windows = []

for line in Path(args.regions).read_text().splitlines():
    chrom, start, end, name = line.split()
    if chrom != "chr20":
        raise SystemExit("Beklenmeyen kromozom.")
    windows.append((
        max(0, int(start) - padding),
        min(chromosome_length, int(end) + padding)
    ))

merged = []
for start, end in sorted(windows):
    if merged and start <= merged[-1][1]:
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    else:
        merged.append((start, end))

sources = {
    "giab": Path(args.giab),
    "deepvariant": Path(args.deepvariant),
    "gatk": Path(args.gatk)
}

folder = Path(args.output_dir)
folder.mkdir(parents=True, exist_ok=True)
targets = {name: folder / f"{name}.padded.vcf" for name in sources}
report = folder / "preparation.json"

if any(p.exists() for p in [*targets.values(), report]):
    raise SystemExit("Hedeflerden biri zaten var; degistirilmedi.")
if any(not p.is_file() for p in sources.values()):
    raise SystemExit("Kaynak VCF dosyalarindan biri eksik.")

counts = {}
for name, source in sources.items():
    selected = 0
    with gzip.open(source, "rt") as reader, targets[name].open("x") as writer:
        for line in reader:
            if line.startswith("#"):
                writer.write(line)
                continue

            f = line.split("\t")
            if f[0] != "chr20":
                continue

            start = int(f[1]) - 1
            end = start + len(f[3])
            if any(start < right and end > left for left, right in merged):
                writer.write(line)
                selected += 1

    counts[name] = selected
    print(name, "| Secilen kayit:", selected, "| Dosya:", targets[name])

with report.open("x") as writer:
    json.dump({
        "chromosome": "chr20",
        "padding_bp": padding,
        "windows_bed_coordinates": merged,
        "sources": {k: str(v) for k, v in sources.items()},
        "record_counts": counts,
        "filtering_applied": False
    }, writer, indent=2)
    writer.write("\n")

print("Hazirlama kaydi:", report)
