#!/usr/bin/env bash
set -euo pipefail

# Run from the project root, regardless of the current directory.
cd "$(dirname "${BASH_SOURCE[0]}")/.."

mode="${1:---check}"
if [[ "$mode" != "--check" && "$mode" != "--run" ]]; then
    echo "Kullanim: bash scripts/01_run_pilot_benchmarks.sh [--check|--run]"
    exit 1
fi

command -v podman >/dev/null
command -v python3 >/dev/null

inputs=(
    docs/software/happy_image_inspect.json
    data/reference/chr20.GRCh38.fa
    data/reference/chr20.GRCh38.fa.fai
    data/benchmark/pilot_chr20_10Mb_11Mb.bed
    data/pilot/giab.pilot.padded.noAD.vcf
    data/pilot/deepvariant.pilot.padded.vcf
    data/pilot/gatk.pilot.padded.noAD.vcf
)

for path in "${inputs[@]}"; do
    if [[ ! -s "$path" ]]; then
        echo "Eksik veya bos dosya: $path"
        exit 1
    fi
done

image_id="$(python3 - <<'PYIMAGE'
import json
from pathlib import Path

path = Path("docs/software/happy_image_inspect.json")
image = json.loads(path.read_text())[0]
assert image.get("Id"), "Imaj kimligi bulunamadi."
print(image["Id"])
PYIMAGE
)"

podman --storage-opt overlay.ignore_chown_errors=true \
    image exists "$image_id"

echo "Girdi dosyalari ve yerel konteyner imaji bulundu."
echo "Imaj: $image_id"
echo "Sinirlar: 1 CPU, konteyner 1 GiB, RTG Java bellegi 512 MiB."

if [[ "$mode" == "--check" ]]; then
    echo "Yalnizca on kontrol yapildi; analiz calistirilmadi."
    exit 0
fi

outdir="results/pilot_workflow"
if [[ -e "$outdir" ]]; then
    echo "Cikti klasoru zaten var; uzerine yazilmadi: $outdir"
    exit 1
fi
mkdir -p "$outdir"

for caller in deepvariant gatk; do
    if [[ "$caller" == "deepvariant" ]]; then
        query="data/pilot/deepvariant.pilot.padded.vcf"
    else
        query="data/pilot/gatk.pilot.padded.noAD.vcf"
    fi

    podman --storage-opt overlay.ignore_chown_errors=true \
        run --rm --cpus 1 --memory 1g \
        -e RTG_MEM=512m \
        -v "$PWD":/work:z -w /work \
        "$image_id" \
        /opt/hap.py/bin/hap.py \
        data/pilot/giab.pilot.padded.noAD.vcf "$query" \
        -r data/reference/chr20.GRCh38.fa \
        -f data/benchmark/pilot_chr20_10Mb_11Mb.bed \
        --engine vcfeval \
        --engine-vcfeval-path /opt/hap.py/libexec/rtg-tools-install/rtg \
        --gender male --threads 1 --pass-only \
        -o "$outdir/$caller" \
        2>&1 | tee "$outdir/$caller.log"
done

echo "Iki pilot karsilastirmasi tamamlandi."
