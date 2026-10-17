#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

inputs=(
    docs/software/happy_image_inspect.json
    results/gene_workflow/deepvariant.vcf.gz
    results/gene_workflow/gatk.vcf.gz
    data/reference/chr20.GRCh38.fa
    data/reference/chr20.GRCh38.fa.fai
    data/regions/chr20_cmrg_evaluable.bed
    data/regions/per_gene/genes.tsv
    data/gene_subset/noAD/giab.vcf
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

image = json.loads(
    Path("docs/software/happy_image_inspect.json").read_text()
)[0]
assert image.get("Id")
print(image["Id"])

tsv = Path("data/regions/per_gene/genes.tsv")
for line in tsv.read_text().splitlines():
    name, filename = line.split("\t")
    bed = tsv.parent / filename
    assert bed.is_file() and bed.stat().st_size > 0, str(bed)
PYIMAGE
)"

podman --storage-opt overlay.ignore_chown_errors=true \
    image exists "$image_id"

outdir="results/gene_workflow/per_gene"
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
        "results/gene_workflow/$caller.vcf.gz" \
        --type ga4gh \
        -r data/reference/chr20.GRCh38.fa \
        -f data/regions/chr20_cmrg_evaluable.bed \
        --adjust-conf-regions data/gene_subset/noAD/giab.vcf \
        --stratification data/regions/per_gene/genes.tsv \
        --threads 1 --no-roc \
        -o "$outdir/$caller" \
        2>&1 | tee "$outdir/$caller.log"
done

echo "Iki program icin gen bazinda olcum tamamlandi."
