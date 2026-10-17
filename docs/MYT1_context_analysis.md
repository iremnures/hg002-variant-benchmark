# MYT1 INDEL error analysis

## Objective

Compare truth-side INDEL recall between existing DeepVariant 1.0.0
and GATK 4.1.8.0 HG002 call sets, and investigate associations with
INDEL length and tandem-repeat overlap.

This is an exploratory analysis of specific historical call sets,
not a general ranking of current variant callers.

## Data and scope

| Component | Definition |
|---|---|
| Sample | HG002 |
| Reference assembly | GRCh38 |
| Benchmark | GIAB HG002 v5.0q small-variant VCF and paired BED |
| Query call sets | DeepVariant 1.0.0; GATK 4.1.8.0 |
| Analysis region | data/regions/per_gene/MYT1.bed |
| Repeat annotation | data/stratifications/GRCh38_AllTandemRepeats.bed.gz |
| Decision source | Annotated hap.py output VCFs in results/gene_workflow/ |
| Analysis unit | GIAB-side VCF INDEL record |

The source CMRG interval spans 139,859 bp. Its intersection with
benchmark regions spans 123,082 bp (88.00%). Performance claims
are restricted to evaluable regions.

The v5.0q benchmark is a draft release. Release documentation is
retained in docs/GIAB_v5.0q_README.md.

## Methods

### Benchmark decisions and record pairing

The analysis reads the TRUTH sample in existing comparison outputs,
selecting BVT=INDEL records with BD=TP or BD=FN.

TP denotes a matched truth record; FN denotes an unmatched truth
record. The script does not rerun hap.py or assign new haplotype
matches. Upstream commands, filtering and engine provenance must
be documented separately.

Cross-output records are paired by chromosome, position, REF and
sorted decoded TRUTH genotype alleles. Decoding allele indices
avoids false mismatches caused by ALT order or unused ALT alleles.
Phase is ignored only for this identity check; existing benchmark
decisions are retained.

Missing genotype components remain explicitly marked and are not
imputed as reference. At chr20:64170884, the original GIAB record
contains a spanning-deletion star allele; both comparison outputs
represent the truth genotype as .|1. This record is retained and
flagged as having a partial truth genotype.

The script stops on duplicate identities, differing paired record
sets, or inconsistent context assignments.

### INDEL length

Length is the maximum absolute length difference between REF and
called literal TRUTH alleles. Bins are 1-5 bp, 6-15 bp and 16+ bp.

Symbolic alleles, spanning-deletion star alleles and missing components
do not contribute to length. A record without a measurable length
change causes the script to stop.

For multiallelic records, the bin describes the largest measurable
called length change. This is not allele-level counting or a
reconstruction of the complete haplotype event.

### Tandem-repeat overlap

Overlap is defined using the record's 0-based, half-open REF span
and the supplied repeat BED. Insertions use the reference anchor.

Any overlap qualifies. This does not require the complete affected
sequence to lie within a repeat, and alternative boundary definitions
have not yet been evaluated.

### Metric

Recall = TRUTH.TP / (TRUTH.TP + TRUTH.FN).

Counts accompany recall. Empty groups have undefined recall, stored
as an empty CSV field and displayed as "No records" in the figure.

FP and precision are outside this truth-side analysis.

## Results

### Paired outcomes

| Outcome | GIAB INDEL records |
|---|---:|
| Matched by both callers | 36 |
| Missed by both callers | 9 |
| Missed only by DeepVariant | 6 |
| Missed only by GATK | 0 |
| Total | 51 |

DeepVariant matched 36/51 records (70.59% recall).
GATK matched 42/51 records (82.35% recall).

### Length and repeat context

Cells show matched / total truth records and recall.

| Length | Repeat overlap | Records | DeepVariant | GATK |
|---|---|---:|---:|---:|
| 1-5 bp | Yes | 13 | 13/13 (100.00%) | 13/13 (100.00%) |
| 1-5 bp | No | 16 | 15/16 (93.75%) | 15/16 (93.75%) |
| 6-15 bp | Yes | 3 | 1/3 (33.33%) | 1/3 (33.33%) |
| 6-15 bp | No | 4 | 2/4 (50.00%) | 4/4 (100.00%) |
| 16+ bp | Yes | 15 | 5/15 (33.33%) | 9/15 (60.00%) |
| 16+ bp | No | 0 | Not estimable | Not estimable |

![MYT1 INDEL recall by length and repeat overlap](../reports/figures/MYT1_indel_recall_by_context.png)

### Interpretation

Both callers recovered all 13 small INDEL records overlapping
tandem repeats. Repeat overlap alone did not imply a missed call.

The lowest observed recall occurred among 16+ bp records, all of
which overlapped repeats. GATK recovered 9/15 and DeepVariant 5/15.

Of the six records missed only by DeepVariant, four were in the
16+ bp repeat-overlapping group and two were in the 6-15 bp
non-overlapping group.

Eight of the nine records missed by both callers overlapped repeats.

## Limitations

- One sample, one selected region and specific historical call sets
  do not establish general caller superiority.
- MYT1 was selected for detailed inspection after reviewing gene-level
  errors. Findings are exploratory, without independent validation.
- Some groups contain only three or four records.
- All 16+ bp records overlap repeats. Independent effects of length
  and repeat context cannot be separated in this dataset.
- Nearby records may share haplotype blocks and should not be treated
  as independent biological events.
- REF-span overlap and maximum allele-length change are operational
  definitions that affect grouping.
- Partial truth genotypes are retained and flagged.
- Raw-read quality, mapping quality and read-level support were not
  assessed. No cause of failure or sample-quality defect is proven.
- TP/FN status does not establish clinical pathogenicity.
- Benchmark coverage and draft-release status limit interpretation.

## Reproduction

Run from the project root after generating gene benchmark outputs:

    python3 scripts/09_analyze_gene_context.py --gene MYT1
    python3 scripts/10_plot_gene_context.py --gene MYT1

Existing outputs are protected against overwrite. Completed outputs
can be inspected without rerunning the scripts.

The analysis script uses the Python standard library. Plotting
requires Matplotlib. The observed working environment was
Python 3.14.7 and Matplotlib 3.10.9.

### Outputs

| Path | Contents |
|---|---|
| results/context_analysis/MYT1/variants.csv | Paired decisions and record-level context |
| results/context_analysis/MYT1/summary.csv | Recall by length and repeat context |
| results/context_analysis/MYT1/metadata.json | Input hashes, script hash and method definitions |
| results/figures/MYT1_indel_recall_by_context.png | Raster figure |
| results/figures/MYT1_indel_recall_by_context.svg | Vector figure |
| results/figures/MYT1_indel_recall_by_context.pdf | PDF figure |

### Validation performed

- The two outputs contained the same 51 paired truth-record identities.
- Counts reproduced the gene-level benchmark: DeepVariant 36 TP / 15 FN;
  GATK 42 TP / 9 FN.
- All 12 saved summary rows reproduced the inspected results.
- Empty groups retained undefined recall instead of zero recall.
- The PNG was visually inspected for labels, counts and empty groups.

## Next validation priorities

1. Review upstream filtering, reference and comparison-engine provenance.
2. Document the AD-removal preprocessing and its rationale.
3. Evaluate sensitivity to partial-genotype records and overlap definitions.
4. Extend selected analyses to other evaluable genes before broader claims.
