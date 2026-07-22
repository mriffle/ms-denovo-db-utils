"""Casanovo mzTab collapsing: sequence normalisation, score handling, ranking."""

from __future__ import annotations

from pathlib import Path

import pytest

from ms_denovo_db_utils.casanovo import (
    MzTabFormatError,
    adjust_score,
    plain_sequence,
    process_files,
)

MZTAB_HEADER = [
    "PSH",
    "PSM_ID",
    "sequence",
    "accession",
    "unique",
    "database",
    "database_version",
    "search_engine",
    "search_engine_score[1]",
    "modifications",
    "retention_time",
    "charge",
    "exp_mass_to_charge",
    "calc_mass_to_charge",
    "spectra_ref",
    "pre",
    "post",
    "start",
    "end",
]


def mztab_file(
    tmp_path: Path,
    rows: list[tuple[str, int, float, float, float]],
    name: str = "run.mztab",
) -> Path:
    """Write an mzTab from (sequence, charge, score, calc_mz, exp_mz) rows."""
    lines = [
        "MTD\tmzTab-version\t1.0.0",
        "MTD\tsoftware[1]\t[MS, MS:1003281, Casanovo, 4.2.1]",
        "\t".join(MZTAB_HEADER),
    ]
    for psm_id, (sequence, charge, score, calc_mz, exp_mz) in enumerate(rows):
        lines.append(
            "\t".join(
                [
                    "PSM",
                    str(psm_id),
                    sequence,
                    "null",
                    "null",
                    "null",
                    "null",
                    "[MS, MS:1003281, Casanovo, 4.2.1]",
                    repr(score),
                    "null",
                    "1200.5",
                    f"{charge}.0",
                    f"{exp_mz:.6f}",
                    f"{calc_mz:.6f}",
                    f"ms_run[1]:index={psm_id}",
                    "null",
                    "null",
                    "null",
                    "null",
                ]
            )
        )
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------
# Sequence normalisation
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("peptidoform", "expected"),
    [
        ("PEPTIDEK", "PEPTIDEK"),
        ("M+15.995DLGEEHFK", "MDLGEEHFK"),
        ("+43.006PEPTIDEK", "PEPTIDEK"),
        ("PEPTIDEK+0.984", "PEPTIDEK"),
        ("-17.027PEPTIDEK", "PEPTIDEK"),
        ("C+57.021PEPTIDEK", "CPEPTIDEK"),
    ],
)
def test_modification_masses_are_stripped(peptidoform: str, expected: str) -> None:
    assert plain_sequence(peptidoform) == expected


def test_stripping_keeps_only_uppercase_residues() -> None:
    assert plain_sequence("pEpTiDeK") == "ETDK"


# --------------------------------------------------------------------------
# Score adjustment
# --------------------------------------------------------------------------
@pytest.mark.parametrize("score", [0.0, 0.5, 0.999, 1.0])
def test_non_negative_scores_are_untouched(score: float) -> None:
    assert adjust_score(score) == score


@pytest.mark.parametrize(("score", "expected"), [(-0.25, 0.75), (-1.0, 0.0), (-0.001, 0.999)])
def test_negative_scores_are_lifted_by_one(score: float, expected: float) -> None:
    """Casanovo marks precursor-filter failures by subtracting 1."""
    assert adjust_score(score) == pytest.approx(expected)


def test_adjusted_scores_land_in_the_unit_interval() -> None:
    assert all(0 <= adjust_score(s) <= 1 for s in (-1.0, -0.5, 0.0, 0.5, 1.0))


# --------------------------------------------------------------------------
# Collapsing PSMs to peptides
# --------------------------------------------------------------------------
def test_highest_scoring_psm_is_kept(tmp_path: Path) -> None:
    path = mztab_file(
        tmp_path,
        [
            ("PEPTIDEK", 2, 0.10, 500.0, 500.0),
            ("PEPTIDEK", 2, 0.95, 500.0, 500.0),
            ("PEPTIDEK", 2, 0.50, 500.0, 500.0),
        ],
    )
    assert process_files([path])["PEPTIDEK"].score == pytest.approx(0.95)


def test_a_negative_score_can_win_after_adjustment(tmp_path: Path) -> None:
    """-0.1 becomes 0.9 and must beat a raw 0.5."""
    path = mztab_file(
        tmp_path, [("PEPTIDEK", 2, 0.5, 500.0, 500.0), ("PEPTIDEK", 2, -0.1, 500.0, 500.0)]
    )
    assert process_files([path])["PEPTIDEK"].score == pytest.approx(0.9)


def test_peptidoforms_of_one_peptide_share_a_row(tmp_path: Path) -> None:
    path = mztab_file(
        tmp_path,
        [
            ("PEPTIDEK", 2, 0.9, 500.0, 500.0),
            ("M+15.995PEPTIDEK", 2, 0.8, 500.0, 500.0),
        ],
    )
    peptides = process_files([path])
    assert set(peptides) == {"PEPTIDEK", "MPEPTIDEK"}


def test_spectra_and_peptidoforms_accumulate_across_files(tmp_path: Path) -> None:
    first = mztab_file(tmp_path, [("PEPTIDEK", 2, 0.5, 500.0, 500.0)], "a.mztab")
    second = mztab_file(tmp_path, [("PEPTIDEK", 3, 0.9, 400.0, 400.0)], "b.mztab")
    peptide = process_files([first, second])["PEPTIDEK"]
    assert peptide.num_spectra == 2
    assert len(peptide.peptidoforms) == 2  # same sequence, different charge
    assert peptide.score == pytest.approx(0.9)


def test_rank_scores_follow_score_order(tmp_path: Path) -> None:
    path = mztab_file(
        tmp_path,
        [
            ("BESTK", 2, 0.99, 500.0, 500.0),
            ("MIDK", 2, 0.50, 500.0, 500.0),
            ("WORSTK", 2, 0.10, 500.0, 500.0),
        ],
    )
    peptides = process_files([path])
    assert peptides["BESTK"].rank_score < peptides["MIDK"].rank_score
    assert peptides["MIDK"].rank_score < peptides["WORSTK"].rank_score


def test_ppm_error_comes_from_the_best_psm(tmp_path: Path) -> None:
    path = mztab_file(
        tmp_path,
        [
            ("PEPTIDEK", 2, 0.10, 500.0, 500.05),
            ("PEPTIDEK", 2, 0.95, 500.0, 500.0),
        ],
    )
    assert process_files([path])["PEPTIDEK"].mz_ppm_error == pytest.approx(0.0)


def test_charge_is_read_from_a_float_field(tmp_path: Path) -> None:
    """Casanovo writes charge as "2.0", which int() alone would reject."""
    path = mztab_file(tmp_path, [("PEPTIDEK", 2, 0.9, 500.0, 500.0)])
    assert process_files([path])["PEPTIDEK"].charge == 2


# --------------------------------------------------------------------------
# Malformed input
# --------------------------------------------------------------------------
def test_a_file_without_a_psh_header_is_an_error(tmp_path: Path) -> None:
    path = tmp_path / "broken.mztab"
    path.write_text("MTD\tmzTab-version\t1.0.0\n")
    with pytest.raises(MzTabFormatError, match="No PSH header"):
        process_files([path])


def test_a_missing_column_is_an_error(tmp_path: Path) -> None:
    path = tmp_path / "broken.mztab"
    path.write_text("PSH\tPSM_ID\tsequence\tcharge\nPSM\t0\tPEPTIDEK\t2\n")
    with pytest.raises(MzTabFormatError, match="Missing expected column"):
        process_files([path])


def test_non_psm_lines_are_ignored(tmp_path: Path) -> None:
    path = mztab_file(tmp_path, [("PEPTIDEK", 2, 0.9, 500.0, 500.0)])
    path.write_text(path.read_text() + "COM\tsome trailing comment\n")
    assert set(process_files([path])) == {"PEPTIDEK"}


def test_a_file_with_no_psms_yields_no_peptides(tmp_path: Path) -> None:
    assert process_files([mztab_file(tmp_path, [])]) == {}
