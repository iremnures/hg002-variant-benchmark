"""Offline tests for input acquisition and reference formatting."""

import hashlib
import io
import json
import runpy
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/24_prepare_inputs.py"
helpers = runpy.run_path(str(SCRIPT), run_name="input_preparation_helpers")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def reference_fixture():
    sequence = b"ACGTACGT"
    expected = b">chr20\nACGT\nACGT\n\n"
    index = b"chr20\t8\t7\t4\t5\n"
    recipe = {
        "reference_accession": "NC_TEST.1",
        "record_name": "chr20",
        "line_bases": 4,
        "newline": "\n",
        "final_newline": True,
        "trailing_blank_lines": 1,
        "sequence_length_bp": 8,
        "expected_sequence_sha256": digest(sequence),
        "expected_fasta_sha256": digest(expected),
        "expected_index_sha256": digest(index),
        "local_fasta_path": "data/reference/chr20.fa",
        "local_index_path": "data/reference/chr20.fa.fai",
    }
    return recipe, expected, index


class InputPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def test_local_source_retrieval_and_verified_install(self):
        source = self.folder / "source.bin"
        source.write_bytes(b"test input\n")
        destination = self.folder / "output/input.bin"

        helpers["install_verified"](
            destination, digest(source.read_bytes()),
            lambda handle: helpers["retrieve"](source.as_uri(), handle),
        )
        self.assertEqual(destination.read_bytes(), source.read_bytes())
        self.assertEqual(list(destination.parent.glob("*.part")), [])

    def test_wrong_checksum_leaves_no_installed_or_partial_file(self):
        destination = self.folder / "output/input.bin"
        with self.assertRaises(ValueError):
            helpers["install_verified"](
                destination, digest(b"expected"),
                lambda handle: handle.write(b"wrong"),
            )
        self.assertFalse(destination.exists())
        self.assertEqual(list(destination.parent.glob("*.part")), [])

    def test_existing_file_is_not_overwritten(self):
        destination = self.folder / "existing.bin"
        destination.write_bytes(b"original")
        with self.assertRaises(FileExistsError):
            helpers["install_verified"](
                destination, digest(b"replacement"),
                lambda handle: handle.write(b"replacement"),
            )
        self.assertEqual(destination.read_bytes(), b"original")
        self.assertEqual(list(self.folder.glob("*.part")), [])

    def test_reference_case_wrapping_header_and_trailing_blank(self):
        recipe, expected, _ = reference_fixture()
        source = io.BytesIO(b">NC_TEST.1 description\nacg\ntacgt\n\n")
        output = io.BytesIO()
        helpers["write_reference"](source, output, recipe)
        self.assertEqual(output.getvalue(), expected)

    def test_wrong_reference_accession_is_rejected(self):
        recipe, _, _ = reference_fixture()
        with self.assertRaises(ValueError):
            helpers["write_reference"](
                io.BytesIO(b">WRONG.1\nACGTACGT\n"),
                io.BytesIO(), recipe,
            )

    def test_wrong_or_truncated_reference_sequence_is_rejected(self):
        recipe, _, _ = reference_fixture()
        for sequence in (b"ACGTTCGT", b"ACGT"):
            with self.subTest(sequence=sequence):
                with self.assertRaises(ValueError):
                    helpers["write_reference"](
                        io.BytesIO(b">NC_TEST.1\n" + sequence + b"\n"),
                        io.BytesIO(), recipe,
                    )

    def make_project(self):
        project = self.folder / "project"
        config = project / "config"
        config.mkdir(parents=True)
        sources = self.folder / "sources"
        sources.mkdir()

        calls = sources / "calls.txt"
        calls.write_bytes(b"tiny callset fixture\n")
        reference_source = sources / "reference.fa"
        reference_source.write_bytes(
            b">NC_TEST.1 source\nacg\ntacgt\n\n"
        )
        recipe, expected_reference, expected_index = reference_fixture()
        manifest = {
            "source_inputs": [
                {
                    "id": "deepvariant",
                    "local_path": "data/calls.txt",
                    "sha256": digest(calls.read_bytes()),
                    "source_url": calls.as_uri(),
                    "source_status": "remote_size_and_md5_verified",
                },
                {
                    "id": "reference_chr20",
                    "local_path": recipe["local_fasta_path"],
                    "sha256": digest(expected_reference),
                    "source_url": reference_source.as_uri(),
                    "source_status": "remote_accession_sequence_sha256_verified",
                },
            ],
            "derived_inputs": [{
                "id": "reference_index",
                "local_path": recipe["local_index_path"],
                "sha256": digest(expected_index),
            }],
        }
        (config / "inputs.json").write_text(json.dumps(manifest))
        (config / "reference_recipe.json").write_text(json.dumps(recipe))
        return project, recipe, expected_reference, expected_index

    def run_setup(self, project, fetch=False):
        command = [sys.executable, str(SCRIPT), "--root", str(project)]
        if fetch:
            command.append("--fetch-missing")
        return subprocess.run(command, capture_output=True, text=True)

    def test_missing_inputs_require_explicit_fetch_option(self):
        project, recipe, _, _ = self.make_project()
        result = self.run_setup(project)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((project / "data/calls.txt").exists())
        self.assertFalse((project / recipe["local_fasta_path"]).exists())

    def test_complete_missing_input_preparation_and_idempotent_recheck(self):
        project, recipe, expected_reference, expected_index = self.make_project()
        first = self.run_setup(project, fetch=True)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)

        self.assertEqual(
            (project / recipe["local_fasta_path"]).read_bytes(),
            expected_reference,
        )
        self.assertEqual(
            (project / recipe["local_index_path"]).read_bytes(),
            expected_index,
        )

        before = {
            str(path.relative_to(project)): path.read_bytes()
            for path in (project / "data").rglob("*") if path.is_file()
        }
        second = self.run_setup(project)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertIn("VERIFIED_EXISTING", second.stdout)

        after = {
            str(path.relative_to(project)): path.read_bytes()
            for path in (project / "data").rglob("*") if path.is_file()
        }
        self.assertEqual(before, after)
        self.assertEqual(list((project / "data").rglob("*.part")), [])

        reports = sorted((project / "results/setup").glob("*.json"))
        self.assertEqual(len(reports), 2)
        self.assertTrue(all(
            json.loads(path.read_text())["status"] == "PASS"
            for path in reports
        ))


if __name__ == "__main__":
    unittest.main()
