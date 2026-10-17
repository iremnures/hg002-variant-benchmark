# Environment setup

## Scope

This document describes environment and benchmark-image setup.
Input acquisition and reference preparation are handled by
scripts/24_prepare_inputs.py. Inputs must be prepared before launching
Nextflow. A fresh-machine environment installation has not yet been tested.

## Observed environment

- Linux with Podman 5.8.4
- Nextflow 26.10.0
- Python 3.14.7
- Matplotlib 3.10.9
- NumPy 2.3.5
- Pillow 11.3.0

The Python requirements record these observed package versions.
They are not a complete transitive dependency lock.

## Python environment

From the project root:

    python3 -m venv .venv
    source .venv/bin/activate
    python -m pip install -r requirements.txt
    python -m unittest discover -s tests -v

Nextflow's native Python processes use python3 from the active PATH.
Keep the virtual environment active when launching the workflow.
The benchmark and quantification commands run inside Podman.

## Benchmark image

The immutable image reference and expected local image ID are recorded
in config/benchmark_container.json. The original tag is informational.

The tag says v0.3.9, while the executable reports
v0.3.8-17-gf15de4a. Both are retained rather than treated as interchangeable.

Use scripts/23_prepare_benchmark_image.py to inspect the image locally
and, when requested, pull it by its recorded digest.

The workflow reads docs/software/happy_image_inspect.json.
The setup script checks compatibility with this record and will not
silently replace a different image identity.

The recorded architecture is included in the container configuration.
Cross-architecture execution has not been validated.

## Resource scope

The workflow schedules one task at a time with one requested CPU.
Benchmark containers have a 1 GiB memory limit and RTG_MEM=512m.

A process memory request in the local Nextflow executor should not be
presented as an independently verified hard memory limit for every
native Python task.

## Prepare the recorded inputs

After activating the Python environment, run from the project root:

    python3 scripts/23_prepare_benchmark_image.py --pull-if-missing
    python3 scripts/24_prepare_inputs.py --fetch-missing

The first command checks or retrieves the digest-pinned benchmark image.
The second verifies existing inputs and retrieves missing ones from
recorded sources. It also prepares the reference FASTA and index.

Existing inputs with unexpected hashes cause failure; they are not
silently replaced. Downloads are prepared in temporary files and
installed only after the expected SHA-256 matches.

The reference recipe preserves the chr20 header, uppercase sequence,
70-base wrapping and one trailing blank line. These details reproduce
the verified reference bytes.

To verify already prepared inputs without downloading:

    python3 scripts/24_prepare_inputs.py

Setup reports are written under results/setup/.

## Test coverage and limits

The sixteen offline tests cover region/context handling and input
preparation, including retrieval from local fixture URLs, checksum
rejection, overwrite protection, reference formatting and repeat execution.

Local fixtures test acquisition logic without network access. They do
not establish that every public endpoint is always available.
