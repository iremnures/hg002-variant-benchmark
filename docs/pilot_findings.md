# HG002 small-variant benchmarking: pilot findings

## Scope and methods

- Sample: HG002; reference coordinates: GRCh38 chromosome 20.
- Pilot interval: chr20:10,000,000–11,000,000 (BED coordinates).
- Evaluated benchmark territory: 988,461 bases.
- Input VCFs include 10 kb of flanking sequence on each side.
- Benchmark: GIAB v5.0q draft small-variant VCF and matching BED.
- Query callsets: DeepVariant 1.0.0 and GATK 4.1.8.0, labelled NovaSeq PCR-free 30x.
- Comparison: hap.py v0.3.8-17-gf15de4a with RTG Tools 3.8.1.
- Contexts: GIAB genome stratifications v3.6, GRCh38 AllTandemRepeats.
- Available FILTER states differ: DeepVariant PASS; GATK '.'.
- In this hap.py version, --pass-only accepts PASS and '.'.

## Results

| Caller | Context | Type | Truth total | FN | FP | Recall (%) | Precision (%) |
|---|---|---|---:|---:|---:|---:|---:|
| DeepVariant | Outside tandem repeats | INDEL | 261 | 14 | 8 | 94.64 | 97.05 |
| DeepVariant | Within tandem repeats | INDEL | 75 | 17 | 13 | 77.33 | 84.15 |
| DeepVariant | Outside tandem repeats | SNP | 1370 | 3 | 2 | 99.78 | 99.85 |
| DeepVariant | Within tandem repeats | SNP | 49 | 5 | 3 | 89.80 | 93.18 |
| GATK | Outside tandem repeats | INDEL | 261 | 29 | 10 | 88.89 | 96.08 |
| GATK | Within tandem repeats | INDEL | 75 | 13 | 16 | 82.67 | 82.22 |
| GATK | Outside tandem repeats | SNP | 1370 | 4 | 5 | 99.71 | 99.63 |
| GATK | Within tandem repeats | SNP | 49 | 6 | 1 | 87.76 | 98.04 |

## Interpretation

Both callsets showed lower SNP and INDEL recall within tandem repeats than outside tandem repeats in this pilot. Within tandem repeats, GATK showed higher INDEL recall, while DeepVariant showed higher INDEL precision.

## Technical checks and limitations

- Removing FORMAT/AD from separate input copies eliminated AD warnings without changing the checked aggregate counts or performance metrics. Original VCFs were retained.
- Overlapping truth records remain; complete representation compatibility has not been established.
- Region-spanning variants can affect more than one context. Context counts must not be summed as independent partitions; the exact INDEL count discrepancy remains to be fully reconciled.
- These are supplied callset comparisons, not a controlled comparison with harmonized filtering.
- The pilot is small and does not establish genome-wide or clinical performance.
- Context association does not establish the cause of individual errors.
- GIAB v5.0q is a draft benchmark, not an error-free truth source.
