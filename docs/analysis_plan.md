# HG002 small-variant benchmarking: draft analysis plan

## Objective
Evaluate small-variant detection performance across selected
genes and genomic contexts in HG002.

## Candidate data
- Benchmark: GIAB HG002 GRCh38 v5.0q small-variant VCF
  and its paired benchmark-region BED.
- Candidate query callsets: published HG002 NovaSeq PCR-free
  30x results from DeepVariant 1.0.0 and GATK 4.1.8.0.
- Genomic context annotations: GIAB v3.6 stratifications.

## Comparison rules
- Use haplotype-aware comparison, following GIAB recommendations.
- Keep non-assessed calls separate from false positives.
- Report SNVs and indels separately.
- Report counts alongside precision, recall and F1.
- Report benchmark coverage within each selected region.
- Select analysis regions before inspecting performance results.

## Checks required before analysis
- Verify exact reference-sequence compatibility.
- Verify input-data provenance and filtering differences.
- Assess DeepVariant training-data overlap.
- Select regions and confirm sufficient evaluable variants.
- Measure runtime and peak memory in a small pilot.

## Interpretation limits
- Results concern one sample and specific historical software versions.
- GIAB v5.0q is an assembly-based draft benchmark.
- Genomic-context associations do not establish error causes.
- This study is not clinical validation.

## Prespecified chr20 gene analysis

- Scope: all seven chr20 genes in the pinned NIST CMRG coordinate file:
  SEMG1, PLTP, CHRNA4, EEF1A2, PTK6, OPRL1 and MYT1.
- This is a limited gene-focused analysis, not comprehensive clinical validation.
- Gene selection uses published coordinates and chromosome availability,
  not observed caller performance.
- All seven genes remain in the coverage report, including unevaluable genes.
- Evaluate whole gene spans intersected with the v5.0q small-variant
  benchmark BED: 188,823 bases in total.
- Include 10 kb flanks in comparison inputs; flanks do not expand the
  intended scoring territory.
- Preserve separate original inputs; remove FORMAT/AD only in working copies.
- Selected REF alleles have been checked against the chr20 reference.
- Report counts and SNP/INDEL metrics separately.
- Undefined metrics must remain undefined; zero truth variants do not
  imply perfect recall.
- Small variant counts limit gene-level comparisons.
- Available FILTER states differ between callsets; results compare
  supplied outputs, not pipelines with harmonized filtering.
- The earlier 1 Mb pilot is exploratory and separate from this gene analysis.
