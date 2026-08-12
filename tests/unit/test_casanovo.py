"""Casanovo mzTab collapsing: sequence validation, score handling, ranking."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pytest

from ms_denovo_db_utils.casanovo import (
    MzTabFormatError,
    adjust_score,
    normalise_modifications,
    peptidoform_key,
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


def mztab_file_with_mods(
    tmp_path: Path,
    rows: Sequence[tuple[str, int, float, float, float, str]],
    name: str = "run.mztab",
) -> Path:
    """Write an mzTab from (sequence, charge, score, calc_mz, exp_mz, mods) rows.

    ``mods`` is the mzTab ``modifications`` field verbatim, so ``"null"`` for an
    unmodified PSM.
    """
    lines = [
        "MTD\tmzTab-version\t1.0.0",
        "MTD\tsoftware[1]\t[MS, MS:1003281, Casanovo, 5.2.0]",
        "\t".join(MZTAB_HEADER),
    ]
    for psm_id, (sequence, charge, score, calc_mz, exp_mz, mods) in enumerate(rows):
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
                    "[MS, MS:1003281, Casanovo, 5.2.0]",
                    repr(score),
                    mods,
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


def mztab_file(
    tmp_path: Path,
    rows: Sequence[tuple[str, int, float, float, float]],
    name: str = "run.mztab",
) -> Path:
    """The same, for PSMs that carry no modification."""
    return mztab_file_with_mods(tmp_path, [(*row, "null") for row in rows], name)


OXIDATION = "1-Oxidation (M):UNIMOD:35"
DEAMIDATION = "4-Deamidated (N):UNIMOD:7"


# --------------------------------------------------------------------------
# Sequence validation -- 5.x writes bare residues, and nothing else is accepted
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "sequence",
    [
        "M+15.995DLGEEHFK",
        "+43.006PEPTIDEK",
        "PEPTIDEK+0.984",
        "-17.027PEPTIDEK",
        "pEpTiDeK",
    ],
)
def test_a_sequence_that_is_not_bare_residues_is_rejected(tmp_path: Path, sequence: str) -> None:
    """Pre-5.x output used to be silently stripped, which is what hid the defect."""
    path = mztab_file(tmp_path, [(sequence, 2, 0.9, 500.0, 500.0)])
    with pytest.raises(MzTabFormatError, match="non-residue character"):
        process_files([path])


def test_the_rejection_names_the_file_and_the_sequence(tmp_path: Path) -> None:
    path = mztab_file(tmp_path, [("M+15.995DLGEEHFK", 2, 0.9, 500.0, 500.0)], "offender.mztab")
    with pytest.raises(MzTabFormatError) as excinfo:
        process_files([path])
    assert "M+15.995DLGEEHFK" in str(excinfo.value)
    assert "offender.mztab" in str(excinfo.value)


# --------------------------------------------------------------------------
# Modification normalisation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("absent", ["null", "", "   "])
def test_an_absent_modification_field_means_unmodified(absent: str) -> None:
    assert normalise_modifications(absent) == ""


def test_modification_order_does_not_change_the_canonical_form() -> None:
    """Sorting makes the key independent of the order the writer chose."""
    assert normalise_modifications(f"{OXIDATION}; {DEAMIDATION}") == normalise_modifications(
        f"{DEAMIDATION}; {OXIDATION}"
    )


def test_different_modification_states_stay_distinct() -> None:
    forms = {
        normalise_modifications(m)
        for m in ("null", OXIDATION, DEAMIDATION, f"{OXIDATION}; {DEAMIDATION}")
    }
    assert len(forms) == 4


def test_the_same_modification_at_a_different_position_is_a_different_form() -> None:
    assert normalise_modifications("1-Oxidation (M):UNIMOD:35") != normalise_modifications(
        "8-Oxidation (M):UNIMOD:35"
    )


def test_the_peptidoform_key_separates_all_three_of_its_parts() -> None:
    keys = {
        peptidoform_key("PEPTIDEK", "null", 2),
        peptidoform_key("PEPTIDEK", "null", 3),
        peptidoform_key("PEPTIDEK", OXIDATION, 2),
        peptidoform_key("PEPTIDEM", "null", 2),
    }
    assert len(keys) == 4


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
    """Modification state does not split a peptide -- it is one row, two peptidoforms."""
    path = mztab_file_with_mods(
        tmp_path,
        [
            ("MPEPTIDEK", 2, 0.9, 500.0, 500.0, "null"),
            ("MPEPTIDEK", 2, 0.8, 500.0, 500.0, OXIDATION),
        ],
    )
    peptides = process_files([path])
    assert set(peptides) == {"MPEPTIDEK"}
    assert len(peptides["MPEPTIDEK"].peptidoforms) == 2


def test_modification_state_is_counted_at_one_charge(tmp_path: Path) -> None:
    """The defect this replaces: same sequence and charge, so only mods separate them.

    Before the modifications column was read these four PSMs counted as one
    peptidoform, because Casanovo 5.x reports an identical `sequence` for all of them.
    """
    path = mztab_file_with_mods(
        tmp_path,
        [
            ("MPEPTIDENK", 2, 0.9, 500.0, 500.0, "null"),
            ("MPEPTIDENK", 2, 0.8, 500.0, 500.0, OXIDATION),
            ("MPEPTIDENK", 2, 0.7, 500.0, 500.0, DEAMIDATION),
            ("MPEPTIDENK", 2, 0.6, 500.0, 500.0, f"{OXIDATION}; {DEAMIDATION}"),
        ],
    )
    peptide = process_files([path])["MPEPTIDENK"]
    assert peptide.num_spectra == 4
    assert len(peptide.peptidoforms) == 4


def test_the_same_peptidoform_seen_twice_counts_once(tmp_path: Path) -> None:
    path = mztab_file_with_mods(
        tmp_path,
        [
            ("MPEPTIDEK", 2, 0.9, 500.0, 500.0, OXIDATION),
            ("MPEPTIDEK", 2, 0.8, 500.0, 500.0, OXIDATION),
        ],
    )
    peptide = process_files([path])["MPEPTIDEK"]
    assert peptide.num_spectra == 2
    assert len(peptide.peptidoforms) == 1


def test_a_reordered_modification_field_is_not_a_second_peptidoform(tmp_path: Path) -> None:
    path = mztab_file_with_mods(
        tmp_path,
        [
            ("MPEPTIDENK", 2, 0.9, 500.0, 500.0, f"{OXIDATION}; {DEAMIDATION}"),
            ("MPEPTIDENK", 2, 0.8, 500.0, 500.0, f"{DEAMIDATION}; {OXIDATION}"),
        ],
    )
    assert len(process_files([path])["MPEPTIDENK"].peptidoforms) == 1


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


def test_a_missing_modifications_column_is_an_error(tmp_path: Path) -> None:
    """Without it every peptidoform count would silently be a charge-state count."""
    header = [c for c in MZTAB_HEADER if c != "modifications"]
    path = tmp_path / "no_mods.mztab"
    path.write_text("\t".join(header) + "\n")
    with pytest.raises(MzTabFormatError, match="Missing expected column"):
        process_files([path])


def test_non_psm_lines_are_ignored(tmp_path: Path) -> None:
    path = mztab_file(tmp_path, [("PEPTIDEK", 2, 0.9, 500.0, 500.0)])
    path.write_text(path.read_text() + "COM\tsome trailing comment\n")
    assert set(process_files([path])) == {"PEPTIDEK"}


def test_a_file_with_no_psms_yields_no_peptides(tmp_path: Path) -> None:
    assert process_files([mztab_file(tmp_path, [])]) == {}
