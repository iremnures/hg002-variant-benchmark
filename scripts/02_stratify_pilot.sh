#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

mode="${1:---check}"
if [[ "$mode" != "--check" && "$mode" != "--run" ]]; then
    echo "Kullanim: bash scripts/02_stratify_pilot.sh [--check|--run]"
    exit 1
fi

command -v podman >/dev/null
command -v python3 >/dev/null

inputs=(
    docs/software/happy_image_inspect.json
    results/pilot_workflow/reproduction_check.json
    results/pilot_workflow/deepvariant.vcf.gz
    results/pilot_workflow/gatk.vcf.gz
    data/reference/chr20.GRCh38.fa
    data/reference/chr20.GRCh38.fa.fai
    data/benchmark/pilot_chr20_10Mb_11Mb.bed
    data/pilot/giab.pilot.padded.noAD.vcf
    data/stratifications/pilot_stratifications.tsv
    data/stratifications/pilot_tandem_repeats.bed
    data/stratifications/pilot_outside_tandem_repeats.bed
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

checks = json.loads(
    Path("results/pilot_workflow/reproduction_check.json").read_text()
)
expected = {
    ("DeepVariant", "SNP"), ("DeepVariant", "INDEL"),
    ("GATK", "SNP"), ("GATK", "INDEL")
}
assert len(checks) == 4
assert {(c["caller"], c["variant_type"]) for c in checks} == expected
assert all(c["passed"] is True for c in checks)

image = json.loads(
    Path("docs/software/happy_image_inspect.json").read_text()
)[0]
assert image.get("Id")
print(image["Id"])
PYIMAGE
)"

podman --storage-opt overlay.ignore_chown_errors=true \
    image exists "$image_id"

if [[ "$mode" == "--check" ]]; then
    echo "Girdiler, onceki kontrol kaydi ve yerel imaj bulundu."
    echo "Bolgesel olcum henuz calistirilmadi."
    exit 0
fi

outdir="results/pilot_workflow/stratified"
if [[ -e "$outdir" ]]; then
    echo "Cikti klasoru zaten var; uzerine yazilmadi: $outdir"
    exit 1
fi
mkdir -p "$outdir"

for caller in deepvariant gatk; do
    podman --storage-opt overlay.ignore_chown_errors=true \
        run --rm --cpus 1 --memory 1g \
        -v "$PWD":/work:z -w /work \
        "$image_id" \
        /opt/hap.py/bin/qfy.py \
        "results/pilot_workflow/$caller.vcf.gz" \
        --type ga4gh \
        -r data/reference/chr20.GRCh38.fa \
        -f data/benchmark/pilot_chr20_10Mb_11Mb.bed \
        --adjust-conf-regions data/pilot/giab.pilot.padded.noAD.vcf \
        --stratification data/stratifications/pilot_stratifications.tsv \
        --threads 1 --no-roc \
        -o "$outdir/$caller" \
        2>&1 | tee "$outdir/$caller.log"
done

echo "Iki program icin bolgesel olcum tamamlandi."
