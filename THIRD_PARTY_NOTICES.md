# Third-party sources and licensing scope

## Original project material

The repository's MIT license covers original project code, workflow
configuration, tests and original documentation. It does not relicense
third-party material, source datasets or external software.

## Retained third-party material

- `docs/GIAB_v5.0q_README.md` is retained upstream GIAB documentation.
  Its attribution and usage statements remain applicable.
- `docs/vcf_headers/` contains metadata extracted from the supplied
  DeepVariant and GATK callsets. These are source records, not original
  project documentation.
- Source metadata and quoted diagnostic messages in provenance records
  remain attributable to their respective providers.

## Data sources

The analysis uses GIAB HG002 v5.0q benchmarks, published HG002 DeepVariant
and GATK callsets, NIST CMRG coordinates, GIAB genome stratifications and
an NCBI chromosome 20 reference sequence.

Source URLs, versions and verification evidence are recorded in
`config/inputs.json`, `config/callsets.json` and `reports/provenance/`.
These records document source identity; they do not grant additional
rights to third-party material.

Original genomic input files are not distributed in this repository.
Users obtaining those files should consult the providers' applicable
usage and attribution statements.

## External software

Nextflow, Python, Podman, hap.py, RTG Tools and Python dependencies are
external software with their own licenses. Their software distributions
and container images are not included in this repository.

The benchmark image identity and observed executable versions are recorded
in `config/benchmark_container.json` and `docs/software/`.

## Derived reports

Tables and figures report this project's analysis of the identified
datasets. Reuse should retain source attribution and the sample, region,
software-version and benchmark limitations documented in the README.
The project license does not replace terms applicable to underlying data.
