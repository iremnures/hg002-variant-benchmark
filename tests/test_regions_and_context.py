"""Behavioral checks for region generation and truth-allele handling."""

import csv
import gzip
import json
import runpy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGION_SCRIPT = ROOT / "scripts/21_prepare_gene_regions.py"
CONTEXT_SCRIPT = ROOT / "scripts/09_analyze_gene_context.py"

regions = runpy.run_path(str(REGION_SCRIPT), run_name="region_test_helpers")
context = runpy.run_path(str(CONTEXT_SCRIPT), run_name="context_test_helpers")


class RegionGenerationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.coordinates = self.folder / "genes.bed"
        self.benchmark = self.folder / "benchmark.bed"
        self.index = self.folder / "reference.fa.fai"
        self.output = self.folder / "output"

        self.coordinates.write_text(
            "chr20\t10\t20\tA\n"
            "chr20\t18\t30\tB\n"
            "chr20\t40\t50\tC\n"
            "chr1\t10\t20\tOTHER\n"
        )
        self.benchmark.write_text(
            "chr20\t0\t10\n"
            "chr20\t12\t15\n"
            "chr20\t15\t18\n"
            "chr20\t19\t25\n"
        )
        self.index.write_text("chr20\t1000\t7\t60\t61\n")

    def run_generator(self):
        return subprocess.run(
            [
                sys.executable, str(REGION_SCRIPT),
                "--coordinates", str(self.coordinates),
                "--benchmark-bed", str(self.benchmark),
                "--reference-index", str(self.index),
                "--output-dir", str(self.output),
            ],
            capture_output=True, text=True,
        )

    def test_half_open_intersection_union_and_zero_coverage(self):
        result = self.run_generator()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

        # [0,10) only touches gene A at 10: it must not be included.
        self.assertEqual(
            (self.output / "per_gene/A.bed").read_text(),
            "chr20\t12\t18\nchr20\t19\t20\n",
        )
        self.assertEqual(
            (self.output / "per_gene/B.bed").read_text(),
            "chr20\t19\t25\n",
        )
        # Shared territory must not be double-counted in the union.
        report = json.loads(
            (self.output / "region_preparation.json").read_text()
        )
        self.assertEqual(report["evaluable_union_bp"], 12)
        self.assertEqual(report["genes_selected"], 3)
        self.assertEqual(report["genes_with_evaluable_regions"], 2)

        with (self.output / "coverage.csv").open(newline="") as handle:
            coverage = {
                row["gene"]: row for row in csv.DictReader(handle)
            }
        self.assertEqual(int(coverage["A"]["evaluable_bp"]), 7)
        self.assertEqual(int(coverage["B"]["evaluable_bp"]), 6)
        self.assertEqual(int(coverage["C"]["evaluable_bp"]), 0)
        self.assertFalse((self.output / "per_gene/C.bed").exists())
        self.assertEqual(
            (self.output / "per_gene/genes.tsv").read_text(),
            "A\tA.bed\nB\tB.bed\n",
        )

    def test_duplicate_gene_is_rejected(self):
        self.coordinates.write_text(
            "chr20\t10\t20\tA\nchr20\t30\t40\tA\n"
        )
        result = self.run_generator()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())

    def test_out_of_bounds_gene_is_rejected(self):
        self.coordinates.write_text("chr20\t990\t1010\tA\n")
        result = self.run_generator()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())

    def test_existing_output_is_preserved(self):
        first = self.run_generator()
        self.assertEqual(first.returncode, 0, first.stderr)
        before = {
            str(path.relative_to(self.output)): path.read_bytes()
            for path in self.output.rglob("*") if path.is_file()
        }

        second = self.run_generator()
        self.assertNotEqual(second.returncode, 0)
        after = {
            str(path.relative_to(self.output)): path.read_bytes()
            for path in self.output.rglob("*") if path.is_file()
        }
        self.assertEqual(before, after)


class ContextAlleleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def write_vcf(self, filename, ref, alt, gt, decision="TP"):
        path = self.folder / filename
        with gzip.open(path, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n"
            )
            handle.write(
                f"chr20\t20\t.\t{ref}\t{alt}\t.\t.\t.\t"
                f"GT:BD:BVT\t{gt}:{decision}:INDEL\n"
            )
        return path

    def test_alt_order_does_not_change_truth_identity(self):
        first = self.write_vcf("first.vcf.gz", "C", "CT,T,CCT", "3|2")
        second = self.write_vcf("second.vcf.gz", "C", "CT,CCT,T", "2|3")

        a = context["read_truth"](first, [("chr20", 0, 100)], [])
        b = context["read_truth"](second, [("chr20", 0, 100)], [])

        self.assertEqual(set(a), set(b))
        self.assertEqual(len(a), 1)
        self.assertEqual(next(iter(a.values()))["max_length"], 2)

    def test_partial_truth_genotype_is_retained_not_imputed(self):
        path = self.write_vcf(
            "partial.vcf.gz", "T", "TATAT", ".|1", decision="FN"
        )
        records = context["read_truth"](path, [("chr20", 0, 100)], [])
        self.assertEqual(len(records), 1)

        key, record = next(iter(records.items()))
        self.assertIn("<MISSING>", key[3])
        self.assertNotIn("T", key[3])
        self.assertTrue(record["partial_gt"])
        self.assertEqual(record["decision"], "FN")
        self.assertEqual(record["max_length"], 4)


class IntervalTests(unittest.TestCase):
    def test_merge_handles_unsorted_nested_and_touching_intervals(self):
        self.assertEqual(
            regions["merge"]([(15, 20), (10, 15), (11, 12), (30, 40)]),
            [(10, 20), (30, 40)],
        )

    def test_half_open_overlap_excludes_boundary_touch(self):
        overlap = context["overlaps"]
        intervals = [("chr20", 10, 20)]
        self.assertFalse(overlap("chr20", 0, 10, intervals))
        self.assertFalse(overlap("chr20", 20, 25, intervals))
        self.assertTrue(overlap("chr20", 19, 21, intervals))
        self.assertFalse(overlap("chr1", 10, 20, intervals))


if __name__ == "__main__":
    unittest.main()
