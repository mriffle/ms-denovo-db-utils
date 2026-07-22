"""End-to-end golden tests over the whole pipeline chain.

The golden files were captured from the pre-refactoring flat scripts, so these
tests are what demonstrate the move into src/ preserved behaviour. They also
pin every later change: any diff here is a change in scientific output and has
to be justified, not just accepted.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ms_denovo_db_utils.cli import (
    build_reset_input,
    collate_into_fasta,
    generate_reverse_decoys,
    process_casanovo_results,
    process_comet_results,
)
from tests.canonical import (
    basename_column,
    canonical_fasta,
    canonical_missing_warning,
    canonical_reset_input,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
GOLDEN = FIXTURES / "golden"

COMET_FILE_COLUMN = 4
CASANOVO_FILE_COLUMN = 3

COMET_RUNS = [str(FIXTURES / "comet/run1.txt"), str(FIXTURES / "comet/run2.txt")]
CASANOVO_RUNS = [str(FIXTURES / "casanovo/run1.mztab"), str(FIXTURES / "casanovo/run2.mztab")]

LIBRARY_DECOY_PREFIX = "LIBRARY_DECOY_"
COMET_DECOY_PREFIX = "COMET_DECOY_"


@pytest.fixture
def comet_peptides(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> Path:
    assert process_comet_results.main(["--decoy_prefix", COMET_DECOY_PREFIX, *COMET_RUNS]) == 0
    path = tmp_path / "comet_peptides.txt"
    path.write_text(capsys.readouterr().out)
    return path


@pytest.fixture
def casanovo_peptides(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> Path:
    assert process_casanovo_results.main(CASANOVO_RUNS) == 0
    path = tmp_path / "casanovo_peptides.txt"
    path.write_text(capsys.readouterr().out)
    return path


def test_comet_peptides_match_golden(comet_peptides: Path) -> None:
    actual = basename_column(comet_peptides.read_text(), COMET_FILE_COLUMN)
    assert actual == (GOLDEN / "comet_peptides.txt").read_text()


def test_casanovo_peptides_match_golden(casanovo_peptides: Path) -> None:
    actual = basename_column(casanovo_peptides.read_text(), CASANOVO_FILE_COLUMN)
    assert actual == (GOLDEN / "casanovo_peptides.txt").read_text()


def test_decoy_generation_matches_golden(capsys: pytest.CaptureFixture[str]) -> None:
    argv = [str(FIXTURES / "fasta/library_targets.fasta"), "--decoy_prefix", LIBRARY_DECOY_PREFIX]
    assert generate_reverse_decoys.main(argv) == 0
    assert capsys.readouterr().out == (GOLDEN / "library_plusdecoys.fasta").read_text()


def test_query_fasta_matches_golden(
    comet_peptides: Path,
    casanovo_peptides: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert collate_into_fasta.main([str(comet_peptides), str(casanovo_peptides)]) == 0
    actual = canonical_fasta(capsys.readouterr().out)
    assert actual == (GOLDEN / "combined_results.fasta").read_text()


def test_reset_input_matches_golden(
    comet_peptides: Path,
    casanovo_peptides: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    argv = [
        str(comet_peptides),
        str(casanovo_peptides),
        str(FIXTURES / "diamond/hits.dmnd.txt"),
        str(GOLDEN / "library_plusdecoys.fasta"),
        LIBRARY_DECOY_PREFIX,
    ]
    assert build_reset_input.main(argv) == 0
    captured = capsys.readouterr()

    assert canonical_reset_input(captured.out) == (GOLDEN / "reset_input.txt").read_text()
    assert canonical_missing_warning(captured.err) == (GOLDEN / "reset_input.stderr").read_text()


def test_reset_input_still_accepts_the_legacy_positional_argument(
    comet_peptides: Path,
    casanovo_peptides: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The deployed nf module passes a comet decoy prefix this tool ignores.

    Dropping the argument outright would break the pipeline against an image
    it has not been updated for, so it stays accepted and unused.
    """
    argv = [
        str(comet_peptides),
        str(casanovo_peptides),
        str(FIXTURES / "diamond/hits.dmnd.txt"),
        str(GOLDEN / "library_plusdecoys.fasta"),
        LIBRARY_DECOY_PREFIX,
        COMET_DECOY_PREFIX,
    ]
    assert build_reset_input.main(argv) == 0
    actual = canonical_reset_input(capsys.readouterr().out)
    assert actual == (GOLDEN / "reset_input.txt").read_text()
