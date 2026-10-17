# Workflow reproducibility validation

## Scope

The Nextflow workflow was executed in a fresh work directory without
`-resume`, starting from existing verified inputs and a local
benchmark container image. All nine tasks completed successfully.

This validates an end-to-end analysis run, not a clean installation
or an automated download-and-rebuild of every input.

## Execution

Console-reported Nextflow version: 26.10.0.

Command, run from the project root:

    NXF_OPTS="-Xms128m -Xmx512m" nextflow run main.nf \
      -work-dir work/clean_run \
      --outdir "$PWD/results/nextflow_clean_run"

The output execution directory was created before launching the run.

## Validation outcome

**PASS: 24 checks passed.**

The validation compared:

- Six prepared VCFs against the established baseline.
- Reference-allele validation results.
- Both callers' benchmark and per-gene summary/extended tables.
- Canonical benchmark VCF record contents.
- MYT1 record-level and summary CSV files.
- Figure-file presence and PNG pixel equality.
- Nine successfully completed, non-cached tasks.

Header timestamps and command paths were excluded from canonical
VCF comparisons. CSV metric comparison uses an absolute tolerance
of 1e-9; other numerical fields are compared exactly.

SVG and PDF file presence was checked; their binary equality was
not required. PNG pixel equality was checked independently of metadata.

## Task measurements

| Task | Task elapsed time | Nextflow-reported peak RSS |
|---|---:|---:|
| PREPARE_GENE_VCFS (chr20 gene subsets) | 14.3s | 20.4 MB |
| CHECK_REF_ALLELES (validate REF alleles) | 77ms | 16.7 MB |
| PREPARE_COMPARISON_INPUTS (prepare AD-free copies) | 96ms | 3.7 MB |
| BENCHMARK_CALLSET (deepvariant) | 7.1s | 55.3 MB |
| BENCHMARK_CALLSET (gatk) | 7.8s | 55.6 MB |
| QUANTIFY_PER_GENE (deepvariant) | 662ms | 3.6 MB |
| ANALYZE_GENE_CONTEXT (MYT1) | 2s | 23.8 MB |
| QUANTIFY_PER_GENE (gatk) | 626ms | 3.6 MB |
| PLOT_GENE_CONTEXT (MYT1) | 1.1s | 117.8 MB |

These are observations from one local run, not performance estimates
across machines or repeated trials. Task times are not total pipeline
wall-clock time. Peak RSS for Podman launcher tasks must not be
presented as verified whole-container or whole-pipeline peak memory.
The workflow benchmarks existing calls; these are not variant-calling
runtime measurements.

## Evidence

- Validator: `scripts/20_validate_complete_workflow.py`
- Machine-readable result: `results/nextflow_clean_run/validation.json`
- Task trace: `results/nextflow_clean_run/execution/trace.tsv`
- Execution report: `results/nextflow_clean_run/execution/report.html`
- Timeline: `results/nextflow_clean_run/execution/timeline.html`

Validated workflow SHA-256:

    dc4f120b4a83d10b5999bd266cdb822928afe27fcdd652f710463adc9dc70236

Validated configuration SHA-256:

    d3b7da3cdd2312aa34ecb95964c124fbd2d10e3120eafb1473acfae9d9c56dbe

## Remaining release work

- Reproduce the derived analysis-region files from source coordinates
  and the benchmark BED.
- Document or automate acquisition and preparation of source inputs.
- Establish portable container retrieval and environment setup.
- Add meaningful automated checks for input and format edge cases.
- Prepare the final README, result overview and repository contents.

Successful reproduction does not independently validate the biological
truth set or remove the limitations described in the analysis reports.

## Final ten-task uncached validation

The final region-generating workflow completed all ten tasks without cache.
All 35 recorded equivalence checks passed.

Source run: results/release_check
Curated evidence: reports/validation/region_workflow.json

The workflow and configuration hashes below were recorded at validation:

    main.nf: 4dab89ca9a28dfe05f91ad494a5e09786ead8b9d7b4c800d832b0ef98e23a89c
    nextflow.config: 26c5dcc6cc628cd23850dcfcc97d463c788c4110e65952c70def81be3adebf93

| Task | Task elapsed time | Nextflow-reported peak RSS |
|---|---:|---:|
| PREPARE_GENE_REGIONS (chr20 CMRG regions) | 105ms | 16.9 MB |
| PREPARE_GENE_VCFS (chr20 gene subsets) | 13.1s | 20.3 MB |
| CHECK_REF_ALLELES (validate REF alleles) | 63ms | 15.9 MB |
| PREPARE_COMPARISON_INPUTS (prepare AD-free copies) | 65ms | 3.8 MB |
| BENCHMARK_CALLSET (deepvariant) | 6.3s | 55.4 MB |
| BENCHMARK_CALLSET (gatk) | 7.4s | 55.4 MB |
| QUANTIFY_PER_GENE (deepvariant) | 552ms | 3.8 MB |
| ANALYZE_GENE_CONTEXT (MYT1) | 1.9s | 23.9 MB |
| QUANTIFY_PER_GENE (gatk) | 565ms | 3.9 MB |
| PLOT_GENE_CONTEXT (MYT1) | 978ms | 16.2 MB |

These measurements are from one local run. Task times are not total
pipeline wall-clock time. Podman launcher RSS is not established
whole-container peak memory.
