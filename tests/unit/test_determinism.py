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

from ms_denovo_db_utils import casanovo, comet

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
    assert comet_num_peptides == "2", "target and decoy Comet hits are both counted"
    assert comet_num_spectra == "2"
    assert "1 decoy peptide(s) on library target rows" in result.stderr


# --------------------------------------------------------------------------
# The other axis: the order the caller passes its input files.
#
# Hash seed is not the only way identical inputs produced different output. Both
# collapsers keep the FIRST best-scoring PSM of a peptide, so a peptide whose best
# score is tied across two files takes its charge, mz_ppm_error and file from
# whichever file was read first. In the pipeline that is the order Nextflow stages
# `.collect()`, which is not fixed between runs -- the drift recorded as A83.
#
# These run in-process rather than through run_cli: the point is argument order,
# and PYTHONHASHSEED has nothing to do with it.
# --------------------------------------------------------------------------
def test_comet_ties_do_not_depend_on_the_order_the_files_are_passed() -> None:
    """Comet caps its e-value at 999, so tied bests are common, not exotic.

    The two fixtures hold the same peptide at the ceiling with different charges and
    different precursor masses, so whichever PSM wins is visible in the output.
    """
    a = FIXTURES / "comet/tied_ceiling_a.txt"
    b = FIXTURES / "comet/tied_ceiling_b.txt"

    forwards = comet.process_files([a, b], COMET_DECOY_PREFIX)
    backwards = comet.process_files([b, a], COMET_DECOY_PREFIX)

    assert list(comet.format_results(forwards)) == list(comet.format_results(backwards))
    # Both PSMs are still counted; only which one represents the peptide was at stake.
    assert forwards["CEILINGPEPTIDEK"].num_spectra == 2


def test_casanovo_ties_do_not_depend_on_the_order_the_files_are_passed() -> None:
    """Latent on the benchmark rather than observed, but a score of exactly 0.0 is
    reachable by underflow, and two files predicting one such peptide would tie."""
    a = FIXTURES / "casanovo/tied_underflow_a.mztab"
    b = FIXTURES / "casanovo/tied_underflow_b.mztab"

    forwards = casanovo.process_files([a, b])
    backwards = casanovo.process_files([b, a])

    assert list(casanovo.format_results(forwards)) == list(casanovo.format_results(backwards))
    assert forwards["UNDERFLOWPEPTLDEK"].score == 0.0
    assert forwards["UNDERFLOWPEPTLDEK"].num_spectra == 2


def test_the_untied_best_psm_still_wins_whatever_the_file_order() -> None:
    """Sorting must not have replaced 'best' with 'first in sorted order'.

    `run2.txt` holds the better SAMPLEPEPTIDER (8e-09 against run1's 1.2e-08) and sorts
    second, so a tool that took the first file's PSM would fail this either way round.
    """
    a = FIXTURES / "comet/run1.txt"
    b = FIXTURES / "comet/run2.txt"

    for order in ([a, b], [b, a]):
        peptides = comet.process_files(order, COMET_DECOY_PREFIX)
        assert peptides["SAMPLEPEPTIDER"].e_value == 8e-09
        assert peptides["SAMPLEPEPTIDER"].source_file == str(b)


@pytest.mark.parametrize("seed", SEEDS)
def test_entrapments_are_identical_under_every_hash_seed(seed: str) -> None:
    """Entrapment generation seeds its shuffles from the peptide sequence, so it must not
    inherit the interpreter's per-process hash salt.

    This is stricter than it looks. The same peptide has to yield the same entrapment in
    the annotated library and in the Comet database, which are generated by separate
    invocations -- so a salted hash would not merely make runs irreproducible, it would
    make the two databases disagree about the same peptide, silently.
    """
    args = [
        str(FIXTURES / "fasta/library_targets.fasta"),
        "--entrapment_prefix",
        "ENTRAPMENT_",
    ]
    first = run_cli("generate_entrapments", args, seed="0")
    assert first.returncode == 0, first.stderr
    result = run_cli("generate_entrapments", args, seed)
    assert result.returncode == 0, result.stderr
    assert result.stdout == first.stdout
