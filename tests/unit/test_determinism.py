"""Identical inputs must produce identical output, run after run.

These use subprocesses because PYTHONHASHSEED only takes effect at interpreter
start. Python randomises string hashing per process, so any code whose result
depends on set or dict-from-set iteration order gives different answers on
different runs -- which for this pipeline means different FDR input.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
GOLDEN = FIXTURES / "golden"

# Seeds chosen to include values that previously produced divergent output and
# aborted runs on the crash fixture.
SEEDS = ["0", "1", "2", "3", "7", "42", "999", "12345"]

LIBRARY_DECOY_PREFIX = "LIBRARY_DECOY_"
COMET_DECOY_PREFIX = "COMET_DECOY_"


def run_cli(module: str, args: list[str], seed: str) -> subprocess.CompletedProcess[str]:
    env = {**os.environ, "PYTHONHASHSEED": seed}
    return subprocess.run(
        [sys.executable, "-m", f"ms_denovo_db_utils.cli.{module}", *args],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


@pytest.fixture(scope="module")
def peptide_tables(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    out = tmp_path_factory.mktemp("tables")
    comet = out / "comet_peptides.txt"
    casanovo = out / "casanovo_peptides.txt"

    result = run_cli(
        "process_comet_results",
        [
            "--decoy_prefix",
            COMET_DECOY_PREFIX,
            str(FIXTURES / "comet/run1.txt"),
            str(FIXTURES / "comet/run2.txt"),
        ],
        seed="0",
    )
    comet.write_text(result.stdout)

    result = run_cli(
        "process_casanovo_results",
        [str(FIXTURES / "casanovo/run1.mztab"), str(FIXTURES / "casanovo/run2.mztab")],
        seed="0",
    )
    casanovo.write_text(result.stdout)
    return comet, casanovo


def build_reset_input_args(comet: Path, casanovo: Path, diamond: str, fasta: str) -> list[str]:
    return [
        str(comet),
        str(casanovo),
        str(FIXTURES / diamond),
        str(fasta),
        LIBRARY_DECOY_PREFIX,
        COMET_DECOY_PREFIX,
    ]


@pytest.mark.parametrize("seed", SEEDS)
def test_reset_input_is_identical_under_every_hash_seed(
    peptide_tables: tuple[Path, Path], seed: str
) -> None:
    comet, casanovo = peptide_tables
    args = build_reset_input_args(
        comet, casanovo, "diamond/hits.dmnd.txt", str(GOLDEN / "library_plusdecoys.fasta")
    )
    result = run_cli("build_reset_input", args, seed)
    assert result.returncode == 0, result.stderr
    assert result.stdout == (GOLDEN / "reset_input.txt").read_text()


@pytest.mark.parametrize("seed", SEEDS)
def test_query_fasta_is_identical_under_every_hash_seed(
    peptide_tables: tuple[Path, Path], seed: str
) -> None:
    """Query order determines DIAMOND's input; it should not drift."""
    comet, casanovo = peptide_tables
    result = run_cli("collate_into_fasta", [str(comet), str(casanovo)], seed)
    assert result.returncode == 0, result.stderr
    assert result.stdout == (GOLDEN / "combined_results.fasta").read_text()


@pytest.mark.parametrize("seed", SEEDS)
def test_missing_peptide_warning_is_identical_under_every_hash_seed(
    peptide_tables: tuple[Path, Path], seed: str
) -> None:
    comet, casanovo = peptide_tables
    args = build_reset_input_args(
        comet, casanovo, "diamond/hits.dmnd.txt", str(GOLDEN / "library_plusdecoys.fasta")
    )
    result = run_cli("build_reset_input", args, seed)
    assert result.stderr == (GOLDEN / "reset_input.stderr").read_text()


@pytest.mark.parametrize("seed", SEEDS)
def test_decoy_comet_hit_on_a_target_group_is_handled_the_same_way_every_run(
    tmp_path_factory: pytest.TempPathFactory, seed: str
) -> None:
    """The fixture that used to abort on some seeds and succeed on others.

    A group holds one Comet decoy hit and one Comet target hit, reached through
    a library decoy and a library target respectively, with the decoy Comet hit
    scoring better. The group's best DIAMOND hit is the target, so the decoy
    Comet hit must always be discarded and never counted.
    """
    out = tmp_path_factory.mktemp(f"crash-{seed}")
    comet = out / "comet_peptides.txt"
    casanovo = out / "casanovo_peptides.txt"

    result = run_cli(
        "process_comet_results",
        ["--decoy_prefix", COMET_DECOY_PREFIX, str(FIXTURES / "comet/crash.txt")],
        seed,
    )
    comet.write_text(result.stdout)
    casanovo.write_text(
        "peptide_sequence\tcharge\tsearch_engine_score[1]\tfile\t"
        "mz_ppm_error\tnum_spectra\trank_score\tnum_peptidoforms\n"
    )

    args = build_reset_input_args(
        comet, casanovo, "diamond/crash.dmnd.txt", str(FIXTURES / "fasta/crash_plusdecoys.fasta")
    )
    result = run_cli("build_reset_input", args, seed)

    assert result.returncode == 0, f"aborted on seed {seed}:\n{result.stderr}"
    rows = [line for line in result.stdout.splitlines()[1:] if line]
    assert len(rows) == 1
    fields = rows[0].split("\t")

    label, comet_num_spectra, comet_num_peptides = fields[1], fields[12], fields[7]
    assert label == "1", "group's best DIAMOND hit is a target"
    assert comet_num_peptides == "1", "the decoy Comet hit must not be counted"
    assert comet_num_spectra == "1"
    assert "Ignoring decoy comet hit for DECOYPEPTIDEK" in result.stderr
