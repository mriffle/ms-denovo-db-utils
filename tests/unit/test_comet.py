"""Comet PSM collapsing: tryptic termini, decoy status, aggregation, ranking."""

from __future__ import annotations

from pathlib import Path

import pytest

from ms_denovo_db_utils.comet import (
    CometFormatError,
    is_c_tryptic,
    is_decoy,
    is_n_tryptic,
    process_files,
)

HEADER = [
    "scan",
    "num",
    "charge",
    "exp_neutral_mass",
    "calc_neutral_mass",
    "e-value",
    "xcorr",
    "delta_cn",
    "sp_score",
    "ions_matched",
    "ions_total",
    "plain_peptide",
    "modified_peptide",
    "prev_aa",
    "next_aa",
    "protein",
    "duplicate_protein_count",
    "modifications",
    "retention_time_sec",
    "sp_rank",
]

DECOY_PREFIX = "COMET_DECOY_"


def comet_file(
    tmp_path: Path,
    rows: list[tuple[str, str, str, float, float, float]],
    name: str = "run.txt",
) -> Path:
    """Write a Comet .txt from (plain, modified, protein, evalue, calc, exp) rows.

    Every row is its own spectrum, ranked 1. Use :func:`ranked_comet_file` to
    write several ranked matches against a single spectrum.
    """
    lines = ["CometVersion 2023.02 rev. 0\tmsms_run_summary", "\t".join(HEADER)]
    for scan, (plain, modified, protein, evalue, calc, exp) in enumerate(rows, start=1):
        lines.append(
            "\t".join(
                [
                    str(scan),
                    "1",
                    "2",
                    f"{exp:.6f}",
                    f"{calc:.6f}",
                    repr(evalue),
                    "3.4",
                    "0.25",
                    "500.1",
                    "12",
                    "24",
                    plain,
                    modified,
                    modified.split(".")[0],
                    modified.split(".")[-1],
                    protein,
                    "0",
                    "-",
                    "1200.5",
                    "1",
                ]
            )
        )
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def ranked_comet_file(
    tmp_path: Path,
    rows: list[tuple[int, int, str, str, float]],
    name: str = "ranked.txt",
    *,
    rank_value: str | None = None,
) -> Path:
    """Write a Comet .txt from explicit (scan, rank, plain, protein, evalue) rows.

    Lets a single spectrum carry several ranked matches, which is what Comet
    emits whenever num_output_lines > 1. `rank_value` overrides the rank column
    verbatim, for testing malformed input.
    """
    lines = ["CometVersion 2023.02 rev. 0\tmsms_run_summary", "\t".join(HEADER)]
    for scan, rank, plain, protein, evalue in rows:
        lines.append(
            "\t".join(
                [
                    str(scan),
                    rank_value if rank_value is not None else str(rank),
                    "2",
                    "1543.741200",
                    "1543.741200",
                    repr(evalue),
                    "3.4",
                    "0.25",
                    "500.1",
                    "12",
                    "24",
                    plain,
                    f"K.{plain}.A",
                    "K",
                    "A",
                    protein,
                    "0",
                    "-",
                    "1200.5",
                    "1",
                ]
            )
        )
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------
# Tryptic termini. Comet's modified_peptide carries flanking residues.
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("modified", "expected"),
    [
        ("K.PEPTIDER.A", 1),
        ("R.PEPTIDER.A", 1),
        ("-.PEPTIDER.A", 1),  # protein N-terminus counts as tryptic
        ("A.PEPTIDER.A", 0),
        ("G.PEPTIDER.A", 0),
        ("M.PEPTIDER.A", 0),
    ],
)
def test_n_terminal_trypticity_reads_the_preceding_residue(modified: str, expected: int) -> None:
    assert is_n_tryptic(modified) == expected


@pytest.mark.parametrize(
    ("plain", "modified", "expected"),
    [
        ("PEPTIDER", "K.PEPTIDER.A", 1),
        ("PEPTIDEK", "K.PEPTIDEK.A", 1),
        ("PEPTIDEA", "K.PEPTIDEA.-", 1),  # protein C-terminus counts as tryptic
        ("PEPTIDEA", "K.PEPTIDEA.G", 0),
        ("PEPTIDEM", "K.PEPTIDEM.G", 0),
    ],
)
def test_c_terminal_trypticity_reads_the_peptide_and_the_next_residue(
    plain: str, modified: str, expected: int
) -> None:
    assert is_c_tryptic(plain, modified) == expected


def test_a_modified_c_terminal_residue_still_counts() -> None:
    """The plain peptide decides, so an inline modification cannot hide R/K."""
    assert is_c_tryptic("PEPTIDER", "K.PEPTIDER[10.0].A") == 1


# --------------------------------------------------------------------------
# Decoy status
# --------------------------------------------------------------------------
def test_single_decoy_protein_is_a_decoy() -> None:
    assert is_decoy(f"{DECOY_PREFIX}sp|X", DECOY_PREFIX) == 1


def test_single_target_protein_is_not_a_decoy() -> None:
    assert is_decoy("sp|X", DECOY_PREFIX) == 0


def test_all_proteins_must_be_decoys_for_the_peptide_to_be_one() -> None:
    assert is_decoy(f"{DECOY_PREFIX}sp|A,{DECOY_PREFIX}sp|B", DECOY_PREFIX) == 1


def test_a_peptide_shared_with_any_target_is_not_a_decoy() -> None:
    """Shared target/decoy peptides must not be scored as decoys."""
    assert is_decoy(f"sp|A,{DECOY_PREFIX}sp|B", DECOY_PREFIX) == 0
    assert is_decoy(f"{DECOY_PREFIX}sp|B,sp|A", DECOY_PREFIX) == 0


def test_decoy_prefix_is_configurable() -> None:
    assert is_decoy("rev_sp|A", "rev_") == 1
    assert is_decoy("rev_sp|A", DECOY_PREFIX) == 0


# --------------------------------------------------------------------------
# Collapsing PSMs to peptides
# --------------------------------------------------------------------------
def test_lowest_evalue_psm_is_kept(tmp_path: Path) -> None:
    path = comet_file(
        tmp_path,
        [
            ("PEPTIDER", "K.PEPTIDER.A", "sp|WORSE", 1e-2, 1000.0, 1000.0),
            ("PEPTIDER", "K.PEPTIDER.A", "sp|BEST", 1e-9, 1000.0, 1000.0),
            ("PEPTIDER", "K.PEPTIDER.A", "sp|MIDDLE", 1e-5, 1000.0, 1000.0),
        ],
    )
    peptides = process_files([path], DECOY_PREFIX)
    assert peptides["PEPTIDER"].e_value == 1e-9
    assert peptides["PEPTIDER"].protein == "sp|BEST"


# --------------------------------------------------------------------------
# Within-spectrum rank. Comet emits num_output_lines matches per spectrum;
# only the best-scoring one is evidence about that spectrum.
# --------------------------------------------------------------------------
def test_runner_up_matches_are_ignored(tmp_path: Path) -> None:
    """A peptide seen only as a spectrum's second-best match is not evidence."""
    path = ranked_comet_file(
        tmp_path,
        [
            (101, 1, "TOPMATCHR", "sp|A", 1e-9),
            (101, 2, "RUNNERUPR", "sp|B", 1e-3),
            (101, 3, "THIRDPLACER", "sp|C", 1e-1),
        ],
    )
    peptides = process_files([path], DECOY_PREFIX)
    assert set(peptides) == {"TOPMATCHR"}


def test_num_spectra_counts_spectra_not_rows(tmp_path: Path) -> None:
    """The bug this guards: with 5 matches per spectrum the count was 5x too high."""
    path = ranked_comet_file(
        tmp_path,
        [
            (101, 1, "PEPTIDER", "sp|A", 1e-9),
            (101, 2, "PEPTIDER", "sp|A", 1e-3),
            (102, 1, "PEPTIDER", "sp|A", 1e-8),
            (102, 2, "PEPTIDER", "sp|A", 1e-2),
        ],
    )
    assert process_files([path], DECOY_PREFIX)["PEPTIDER"].num_spectra == 2


def test_peptidoforms_from_runner_up_matches_are_excluded(tmp_path: Path) -> None:
    path = ranked_comet_file(
        tmp_path,
        [
            (101, 1, "PEPTIDER", "sp|A", 1e-9),
            (101, 2, "PEPTIDEK", "sp|A", 1e-3),
        ],
    )
    peptides = process_files([path], DECOY_PREFIX)
    assert len(peptides["PEPTIDER"].peptidoforms) == 1
    assert "PEPTIDEK" not in peptides


def test_tied_top_ranked_matches_are_both_kept(tmp_path: Path) -> None:
    """Both are genuinely top-scoring; dropping one would be an arbitrary choice."""
    path = ranked_comet_file(
        tmp_path,
        [
            (101, 1, "FIRSTTIEDR", "sp|A", 1e-9),
            (101, 1, "SECONDTIEDR", "sp|B", 1e-9),
        ],
    )
    assert set(process_files([path], DECOY_PREFIX)) == {"FIRSTTIEDR", "SECONDTIEDR"}


def test_a_file_without_a_rank_column_still_parses(tmp_path: Path) -> None:
    """`num` is optional: older or hand-made files keep their previous behaviour."""
    path = ranked_comet_file(tmp_path, [(101, 1, "PEPTIDER", "sp|A", 1e-9)])
    text = path.read_text().splitlines()
    header = text[1].split("\t")
    drop = header.index("num")
    stripped = [text[0]]
    stripped += [
        "\t".join(c for i, c in enumerate(line.split("\t")) if i != drop) for line in text[1:]
    ]
    path.write_text("\n".join(stripped) + "\n")

    assert set(process_files([path], DECOY_PREFIX)) == {"PEPTIDER"}


def test_a_non_numeric_rank_is_an_error(tmp_path: Path) -> None:
    path = ranked_comet_file(tmp_path, [(101, 1, "PEPTIDER", "sp|A", 1e-9)], rank_value="best")
    with pytest.raises(CometFormatError, match="Non-numeric"):
        process_files([path], DECOY_PREFIX)


def test_every_psm_counts_toward_num_spectra(tmp_path: Path) -> None:
    path = comet_file(
        tmp_path,
        [("PEPTIDER", "K.PEPTIDER.A", "sp|A", e, 1000.0, 1000.0) for e in (1e-2, 1e-9, 1e-5)],
    )
    assert process_files([path], DECOY_PREFIX)["PEPTIDER"].num_spectra == 3


def test_spectra_counts_accumulate_across_files(tmp_path: Path) -> None:
    """Nextflow passes every run's Comet output in a single invocation."""
    first = comet_file(
        tmp_path, [("PEPTIDER", "K.PEPTIDER.A", "sp|A", 1e-5, 1000.0, 1000.0)], "run1.txt"
    )
    second = comet_file(
        tmp_path, [("PEPTIDER", "K.PEPTIDER.A", "sp|A", 1e-9, 1000.0, 1000.0)], "run2.txt"
    )
    peptides = process_files([first, second], DECOY_PREFIX)
    assert peptides["PEPTIDER"].num_spectra == 2
    assert peptides["PEPTIDER"].e_value == 1e-9
    assert peptides["PEPTIDER"].source_file == str(second)


def test_peptidoforms_are_counted_distinctly(tmp_path: Path) -> None:
    path = comet_file(
        tmp_path,
        [
            ("PEPTIDER", "K.PEPTIDER.A", "sp|A", 1e-5, 1000.0, 1000.0),
            ("PEPTIDER", "K.PEPTIDER.A", "sp|A", 1e-6, 1000.0, 1000.0),
            ("PEPTIDER", "K.PEPT[79.97]IDER.A", "sp|A", 1e-7, 1000.0, 1000.0),
        ],
    )
    assert len(process_files([path], DECOY_PREFIX)["PEPTIDER"].peptidoforms) == 2


def test_the_same_peptidoform_at_a_different_charge_counts_separately(tmp_path: Path) -> None:
    """Charge is part of the peptidoform identity; commit 87a90d3 made it so."""
    path = comet_file(tmp_path, [("PEPTIDER", "K.PEPTIDER.A", "sp|A", 1e-5, 1000.0, 1000.0)])
    text = path.read_text().replace("\n", "\n", 1).splitlines()
    duplicated = [*text, text[-1].replace("\t2\t", "\t3\t", 1)]
    path.write_text("\n".join(duplicated) + "\n")
    assert len(process_files([path], DECOY_PREFIX)["PEPTIDER"].peptidoforms) == 2


def test_rank_scores_follow_evalue_order(tmp_path: Path) -> None:
    path = comet_file(
        tmp_path,
        [
            ("BESTK", "K.BESTK.A", "sp|A", 1e-9, 1000.0, 1000.0),
            ("MIDK", "K.MIDK.A", "sp|A", 1e-5, 1000.0, 1000.0),
            ("WORSTK", "K.WORSTK.A", "sp|A", 1e-1, 1000.0, 1000.0),
        ],
    )
    peptides = process_files([path], DECOY_PREFIX)
    assert peptides["BESTK"].rank_score < peptides["MIDK"].rank_score
    assert peptides["MIDK"].rank_score < peptides["WORSTK"].rank_score


def test_tied_evalues_receive_the_same_rank(tmp_path: Path) -> None:
    path = comet_file(
        tmp_path,
        [
            ("ONEK", "K.ONEK.A", "sp|A", 5e-6, 1000.0, 1000.0),
            ("TWOK", "K.TWOK.A", "sp|A", 5e-6, 1000.0, 1000.0),
        ],
    )
    peptides = process_files([path], DECOY_PREFIX)
    assert peptides["ONEK"].rank_score == peptides["TWOK"].rank_score


def test_an_evalue_of_zero_is_accepted(tmp_path: Path) -> None:
    """Comet does report 0.0; it must not blow up the transform downstream."""
    path = comet_file(tmp_path, [("PEPTIDER", "K.PEPTIDER.A", "sp|A", 0.0, 1000.0, 1000.0)])
    assert process_files([path], DECOY_PREFIX)["PEPTIDER"].e_value == 0.0


def test_ppm_error_is_derived_from_the_best_psm(tmp_path: Path) -> None:
    path = comet_file(
        tmp_path,
        [
            ("PEPTIDER", "K.PEPTIDER.A", "sp|A", 1e-2, 1000.0, 1000.02),
            ("PEPTIDER", "K.PEPTIDER.A", "sp|A", 1e-9, 1000.0, 1000.0),
        ],
    )
    assert process_files([path], DECOY_PREFIX)["PEPTIDER"].mz_ppm_error == pytest.approx(0.0)


# --------------------------------------------------------------------------
# Malformed input
# --------------------------------------------------------------------------
def test_a_missing_column_is_an_error_not_an_empty_file(tmp_path: Path) -> None:
    """Silently emitting nothing would read as success to Nextflow."""
    path = tmp_path / "broken.txt"
    path.write_text("CometVersion 2023.02\nscan\tcharge\tprotein\n1\t2\tsp|A\n")
    with pytest.raises(CometFormatError, match="Missing expected column"):
        process_files([path], DECOY_PREFIX)


def test_a_file_with_no_psms_yields_no_peptides(tmp_path: Path) -> None:
    path = comet_file(tmp_path, [])
    assert process_files([path], DECOY_PREFIX) == {}
