# DeepVariant training-overlap assessment

## Purpose

Distinguish sample-level training overlap from overlap with the
specific genomic regions evaluated in this project.

## Evidence

The official historical training-data document, retained under the
r1.1 tag, describes v1.0 WGS training using HG002-HG004 examples from
the additional v4.2 confident regions relative to v3.3.2.
It also states that chromosomes 20-22 were excluded from training.

Source:
https://github.com/google/deepvariant/blob/r1.1/docs/deepvariant-details-training-data.md

The official 2021 historical overview likewise describes chromosome 20
as excluded from training for its WGS/PacBio evaluation.

Source:
https://google.github.io/deepvariant/posts/2021-06-08-deepvariant-over-the-years/

These are historical descriptions of released models, not proof of
the checkpoint used to generate the query VCF.

## Project-specific evidence

- The query VCF header identifies DeepVariant version 1.0.0.
- The compressed query VCF matches the public source object's size
  and MD5; see results/validation/callset_source_verification.json.
- Evaluated regions are on chromosome 20.
- The exact checkpoint and training manifest for this query VCF have
  not been established.

## Interpretation

Sample-level independence must not be claimed.

The documented chromosome holdout supports the choice of chr20,
conditional on use of a model following that documented policy.
It does not establish complete independence from model selection,
tuning or development decisions.

Neither the later benchmark release date nor exact file identity
alone proves absence of training overlap.

## Reporting statement

This project evaluates a verified public DeepVariant 1.0.0 call set
on chr20. Historical release documentation describes chr20 as held
out from training, while HG002 contributed training examples elsewhere.
The exact checkpoint used for the call set remains unverified.

## Follow-up

Seek the dataset's generation recipe or model/checkpoint provenance.
If unavailable, retain this uncertainty in the final report.

This assessment does not modify variant calls or benchmark results.
