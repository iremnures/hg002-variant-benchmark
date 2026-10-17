# Gene-region benchmark methods and preprocessing

## Scope

This project evaluates existing HG002 DeepVariant 1.0.0 and
GATK 4.1.8.0 small-variant call sets against the GIAB HG002
GRCh38 v5.0q benchmark and its paired confidence-region BED.

It compares specific historical call sets. It does not establish
general superiority of current caller versions.

## Comparison configuration

The command is retained in scripts/07_run_gene_benchmarks.sh.

| Setting | Configuration |
|---|---|
| Comparison framework | hap.py v0.3.8-17-gf15de4a |
| Comparison engine | RTG Tools 3.8.1 / vcfeval |
| Reference | data/reference/chr20.GRCh38.fa |
| Truth VCF | data/gene_subset/noAD/giab.vcf |
| Query VCFs | data/gene_subset/noAD/{deepvariant,gatk}.vcf |
| Confidence regions | data/regions/chr20_cmrg_evaluable.bed |
| Query filtering option | --pass-only |
| Sample sex setting | --gender male |
| Threads | 1 |
| Container memory limit | 1 GiB |
| RTG Java memory setting | RTG_MEM=512m |

The local container image is selected by its recorded image ID in
docs/software/happy_image_inspect.json. The executed image reports hap.py v0.3.8-17-gf15de4a and
RTG Tools 3.8.1. Raw version-query outputs are retained in
docs/software/benchmark_versions.json.
A portable image retrieval reference remains to be documented.

## Filtering policy

Inspection of the installed hap.py/pre.py code showed that
--pass-only sets the preprocessing filter policy to "*".
The preprocessing code includes both PASS and "." among allowed
FILTER values.

These values have different VCF meanings: PASS reports that filters
were passed, whereas "." reports that filters were not applied.
Acceptance by preprocessing does not make them equivalent evidence
of call quality.

Input FILTER distributions in the padded gene-region subsets were:

| Dataset | FILTER | Records |
|---|---|---:|
| GIAB | . | 1105 |
| GIAB | GAP1 | 8 |
| DeepVariant | PASS | 842 |
| DeepVariant | RefCall | 1354 |
| GATK | . | 801 |

These are padded-input counts, not evaluated TP/FN counts.

DeepVariant RefCall records are excluded under the configured query
filter policy. GATK records have FILTER="."; no additional GATK
quality-filtering step is shown in the benchmark command.

Consequently, this is not a comparison at harmonized quality
thresholds. That limitation must accompany caller comparisons.

The GIAB header defines GAP1 as "Uncalled in the first haplotype".
The run does not enable --usefiltered-truth; named filtered truth
records are excluded under the configured default truth policy.

ALL and PASS summary rows originate from the same pass-only
comparison and are not independent filtered/unfiltered experiments.

## AD compatibility preprocessing

### Observed issue

Historical pilot logs contain AD cardinality warnings, including:

    too many AD fields at chr20:10028312 max_ad = 2 retrieved: 3

Evidence: logs/deepvariant_pilot.log, line 7.

Additional warnings occur in logs/gatk_pilot_noAD.log, including:

    too many AD fields at chr20:10034909 max_ad = 2 retrieved: 3

Evidence: logs/gatk_pilot_noAD.log, line 9.

These are warnings, not proof of a failed run. A filename containing
"noAD" does not establish which inputs were modified.

Previous pilot inspection associated the issue with multiallelic
records containing spanning-deletion star alleles. The exact internal
software mechanism was not independently established.

### Implemented workaround

scripts/06_prepare_gene_benchmark_inputs.py removes FORMAT/AD and
its header definition from comparison copies of all three inputs.

Original padded VCFs are preserved. No quality or genotype filter
is implemented in this preparation script.

Before preparation, the script requires nonzero reference-check
counts and no reported REF mismatches. It also checks that output
record counts match the reference-check report.

| Input | Records retained | Records with AD removed |
|---|---:|---:|
| GIAB | 1113 | 1113 |
| DeepVariant | 2196 | 2196 |
| GATK | 801 | 801 |

The preparation report is stored at:
data/gene_subset/noAD/preparation.json

### Content-preservation validation

A separate source-versus-output comparison verified:

- Headers were unchanged except for removal of the AD definition.
- Record order and counts were preserved.
- The first eight VCF columns were unchanged.
- FORMAT keys were unchanged except for AD removal.
- Remaining sample values, including GT, were preserved.
- Missing trailing FORMAT values were compared semantically as ".".

This establishes preservation of non-AD content. It does not by
itself prove identical benchmark outcomes between the warning-
producing and AD-stripped runs.

The targeted search for "too many AD fields", "max_ad" and
"retrieved:" found matches in historical pilot logs, with no matches
returned from the gene_workflow directory in that inspection.

## Per-gene quantification

scripts/08_quantify_per_gene.sh runs qfy.py on the existing
comparison VCFs; it does not rerun variant calling or vcfeval.

It uses the reference, evaluable BED and:
data/regions/per_gene/genes.tsv

The command also supplies --adjust-conf-regions with the prepared
truth VCF. Its exact region-adjustment behavior should be documented
from the installed version before final methods publication.

## Linked analyses

- [MYT1 context analysis](MYT1_context_analysis.md)

## Remaining reproducibility work

1. Document portable image retrieval provenance.
2. Document source URLs, checksums and query/reference provenance.
3. Review the confidence-region adjustment used during quantification.
4. Assess preprocessing, engine/version and context-definition sensitivity.
5. Add workflow orchestration, dependency records and repository guidance.

## Executable validation and software notes

The non-AD content-preservation check is implemented in
scripts/11_validate_ad_preprocessing.py.

All three inputs passed. The machine-readable report, including file
hashes and checked record counts, is retained in:
results/validation/ad_content_preservation.json

Both container software version queries returned exit code zero.
The reference-location warnings occurred during version-only queries,
without a mounted project reference. The benchmark command separately
supplies the reference through its explicit -r option.

The evaluation software is historical. These results describe the
recorded environment; agreement with a different evaluation version
or engine has not yet been assessed.
