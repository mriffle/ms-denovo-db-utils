"""The image must produce the same output as the source tree.

These run the tools through the exact absolute paths and command forms that
nf-ms-denovo-db uses, and compare against the same golden files the unit tests
assert on. A broken install layer, a missing module or a shim that silently
fails would show up here and nowhere else.
"""

from __future__ import annotations

import subprocess
from collections.abc import Callable
from pathlib import Path

import pytest

from tests.canonical import basename_column

pytestmark = pytest.mark.docker

Runner = Callable[..., "subprocess.CompletedProcess[str]"]

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
GOLDEN = FIXTURES / "golden"

COMET_FILE_COLUMN = 4
CASANOVO_FILE_COLUMN = 3

# The pipeline's own invocation style: bash with process substitution through
# tee, which is what modules/create_peptide_fasta.nf and friends emit.
PIPELINE_CHAIN = """
set -euo pipefail
cd /out
python3 /usr/local/bin/process_comet_results.py --decoy_prefix COMET_DECOY_ \
    /data/comet/run1.txt /data/comet/run2.txt \
    > >(tee comet_peptides.txt) 2> >(tee comet.stderr >&2)
python3 /usr/local/bin/process_casanovo_results.py \
    /data/casanovo/run1.mztab /data/casanovo/run2.mztab \
    > >(tee casanovo_peptides.txt) 2> >(tee casanovo.stderr >&2)
python3 /usr/local/bin/collate_into_fasta.py comet_peptides.txt casanovo_peptides.txt \
    > >(tee combined_results.fasta) 2> >(tee collate.stderr >&2)
python3 /usr/local/bin/generate_reverse_decoys.py /data/fasta/library_targets.fasta \
    --decoy_prefix LIBRARY_DECOY_ > library_plusdecoys.fasta
python3 /usr/local/bin/build_reset_input.py \
    comet_peptides.txt casanovo_peptides.txt /data/diamond/hits.dmnd.txt \
    library_plusdecoys.fasta LIBRARY_DECOY_ COMET_DECOY_ \
    > >(tee reset_input.txt) 2> >(tee reset_input.stderr >&2)
sleep 0.3
"""


@pytest.fixture(scope="module")
def pipeline_output(
    run_in_container: Runner,
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    """Run the whole chain inside the container and return the output dir."""
    out_dir = tmp_path_factory.mktemp("container-out")
    result = run_in_container(
        ["bash", "-c", PIPELINE_CHAIN],
        mounts={FIXTURES: "/data", out_dir: "/out"},
        as_current_user=True,
    )
    if result.returncode != 0:
        pytest.fail(f"pipeline chain failed in container:\n{result.stdout}\n{result.stderr}")
    return out_dir


def test_comet_peptides_match_golden(pipeline_output: Path) -> None:
    actual = basename_column(
        (pipeline_output / "comet_peptides.txt").read_text(), COMET_FILE_COLUMN
    )
    assert actual == (GOLDEN / "comet_peptides.txt").read_text()


def test_casanovo_peptides_match_golden(pipeline_output: Path) -> None:
    actual = basename_column(
        (pipeline_output / "casanovo_peptides.txt").read_text(), CASANOVO_FILE_COLUMN
    )
    assert actual == (GOLDEN / "casanovo_peptides.txt").read_text()


def test_query_fasta_matches_golden(pipeline_output: Path) -> None:
    actual = (pipeline_output / "combined_results.fasta").read_text()
    assert actual == (GOLDEN / "combined_results.fasta").read_text()


def test_decoy_library_matches_golden(pipeline_output: Path) -> None:
    actual = (pipeline_output / "library_plusdecoys.fasta").read_text()
    assert actual == (GOLDEN / "library_plusdecoys.fasta").read_text()


def test_reset_input_matches_golden(pipeline_output: Path) -> None:
    actual = (pipeline_output / "reset_input.txt").read_text()
    assert actual == (GOLDEN / "reset_input.txt").read_text()


def test_reset_input_warning_matches_golden(pipeline_output: Path) -> None:
    actual = (pipeline_output / "reset_input.stderr").read_text()
    assert actual == (GOLDEN / "reset_input.stderr").read_text()


def test_gzipped_library_yields_gzipped_output(run_in_container: Runner, tmp_path: Path) -> None:
    """GENERATE_LIBRARY_DECOYS names its output .gz based on this behaviour."""
    result = run_in_container(
        [
            "bash",
            "-c",
            "python3 /usr/local/bin/generate_reverse_decoys.py "
            "/data/fasta/library_targets.fasta.gz --decoy_prefix LIBRARY_DECOY_ "
            "> /out/out.fasta.gz",
        ],
        mounts={FIXTURES: "/data", tmp_path: "/out"},
        as_current_user=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "out.fasta.gz").read_bytes()[:2] == b"\x1f\x8b"
