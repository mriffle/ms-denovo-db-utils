"""I == L merging, applied where the two engines' results are joined.

Casanovo emits only L; Comet reports whatever the database spells. Without this the
I-form and the L-form of one observed peptide are two hypotheses in the FDR pool.
The merge happens at the join and NOT in the query FASTA, so both spellings still
align against the un-normalised library at their own best spelling.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from ms_denovo_db_utils.cli import build_reset_input
from ms_denovo_db_utils.diamond import add_subject_sequences, read_hits
from ms_denovo_db_utils.peptide import normalise_il
from ms_denovo_db_utils.reset_input import (
    OUTPUT_COLUMNS,
    read_casanovo_peptides,
    read_comet_peptides,
    write_reset_input,
)

DECOY = "LIBRARY_DECOY_"

COMET_HEADER = (
    "plain_peptide\tcharge\te-value\tprotein\tfile\ttryptic_n\ttryptic_c\t"
    "num_spectra\tmz_ppm_error\tis_decoy\tproteins\trank_score\tnum_peptidoforms"
)
CASANOVO_HEADER = (
    "peptide_sequence\tcharge\tsearch_engine_score[1]\tfile\tmz_ppm_error\t"
    "num_spectra\trank_score\tnum_peptidoforms"
)


def comet_table(tmp_path: Path, peptides: list[str]) -> Path:
    path = tmp_path / "comet_peptides.txt"
    rows = [COMET_HEADER] + [
        f"{p}\t2\t1e-09\tsp|Q1|P\trun.txt\t1\t1\t1\t1.00\t0\tsp|Q1|P\t0.5\t1" for p in peptides
    ]
    path.write_text("\n".join(rows) + "\n")
    return path


def casanovo_table(tmp_path: Path, peptides: list[str]) -> Path:
    path = tmp_path / "casanovo_peptides.txt"
    rows = [CASANOVO_HEADER] + [f"{p}\t2\t0.9\trun.mztab\t2.00\t1\t0.5\t1" for p in peptides]
    path.write_text("\n".join(rows) + "\n")
    return path


def diamond_table(tmp_path: Path, rows: list[tuple[str, str, int, int, float]]) -> Path:
    """(query, subject, sstart, send, evalue) -> a DIAMOND --outfmt 6 file."""
    path = tmp_path / "hits.dmnd.txt"
    path.write_text(
        "\n".join(
            f"{q}\t{s}\t100.0\t{send - sstart + 1}\t0\t0\t1\t{send - sstart + 1}\t"
            f"{sstart}\t{send}\t{ev}\t50.0"
            for q, s, sstart, send, ev in rows
        )
        + "\n"
    )
    return path


# --------------------------------------------------------------------------
# The canonical spelling
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("PEPTIDEK", "PEPTLDEK"),
        ("PEPTLDEK", "PEPTLDEK"),
        ("IIIK", "LLLK"),
        ("NOTHINGHERE", "NOTHLNGHERE"),
    ],
)
def test_isoleucine_becomes_leucine(raw: str, expected: str) -> None:
    assert normalise_il(raw) == expected


def test_normalisation_is_idempotent() -> None:
    assert normalise_il(normalise_il("PEPTIDEIK")) == normalise_il("PEPTIDEIK")


# --------------------------------------------------------------------------
# read_hits: the two spellings compete, and ambiguity is judged on their union
# --------------------------------------------------------------------------
def test_the_two_spellings_become_one_query(tmp_path: Path) -> None:
    path = diamond_table(
        tmp_path,
        [("PEPTIDEK", "sp|P1|A", 1, 8, 1e-09), ("PEPTLDEK", "sp|P1|A", 1, 8, 1e-05)],
    )
    assert set(read_hits(path, DECOY, il_equivalence=True)) == {"PEPTLDEK"}
    assert set(read_hits(path, DECOY)) == {"PEPTIDEK", "PEPTLDEK"}


def test_the_better_evalue_wins_across_spellings(tmp_path: Path) -> None:
    path = diamond_table(
        tmp_path,
        [("PEPTIDEK", "sp|P1|A", 1, 8, 1e-09), ("PEPTLDEK", "sp|P2|B", 1, 8, 1e-05)],
    )
    hit = read_hits(path, DECOY, il_equivalence=True)["PEPTLDEK"]
    assert hit.sseqid == "sp|P1|A"
    assert hit.qseqid == "PEPTIDEK", "the surviving hit keeps the spelling that aligned"


def test_ambiguity_is_judged_on_the_merged_peptide(tmp_path: Path) -> None:
    """Neither spelling spans the axis alone; together they do, so the merge is dropped.

    This is why the merge happens as the file is read rather than by collapsing two
    already-filtered maps -- that would let this peptide through.
    """
    path = diamond_table(
        tmp_path,
        [("PEPTIDEK", "sp|P1|A", 1, 8, 1e-09), ("PEPTLDEK", f"{DECOY}sp|P1|A", 1, 8, 1e-08)],
    )
    assert read_hits(path, DECOY, il_equivalence=True) == {}
    assert set(read_hits(path, DECOY)) == {"PEPTIDEK", "PEPTLDEK"}


# --------------------------------------------------------------------------
# The engine tables
# --------------------------------------------------------------------------
def test_engine_tables_are_keyed_on_the_canonical_spelling(tmp_path: Path) -> None:
    comet = read_comet_peptides(comet_table(tmp_path, ["PEPTIDEK"]), il_equivalence=True)
    casa = read_casanovo_peptides(casanovo_table(tmp_path, ["PEPTLDEK"]), il_equivalence=True)
    assert set(comet) == set(casa) == {"PEPTLDEK"}


def test_without_the_flag_the_tables_keep_their_raw_spelling(tmp_path: Path) -> None:
    comet = read_comet_peptides(comet_table(tmp_path, ["PEPTIDEK"]))
    assert set(comet) == {"PEPTIDEK"}


@pytest.mark.parametrize("reader", [read_comet_peptides, read_casanovo_peptides])
def test_two_rows_of_one_table_colliding_is_an_error(tmp_path: Path, reader: object) -> None:
    """Isobaric database peptides can collide; merging their peptidoform counts is not
    possible from these files, so this refuses rather than guessing."""
    table = (
        comet_table(tmp_path, ["PEPTIDEK", "PEPTLDEK"])
        if reader is read_comet_peptides
        else casanovo_table(tmp_path, ["PEPTIDEK", "PEPTLDEK"])
    )
    with pytest.raises(ValueError, match="both normalise to"):
        reader(table, il_equivalence=True)  # type: ignore[operator]


def test_the_collision_error_names_both_peptides(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as excinfo:
        read_comet_peptides(comet_table(tmp_path, ["PEPTIDEK", "PEPTLDEK"]), il_equivalence=True)
    assert "PEPTIDEK" in str(excinfo.value)
    assert "PEPTLDEK" in str(excinfo.value)


def test_a_collision_is_fine_without_the_flag(tmp_path: Path) -> None:
    assert len(read_comet_peptides(comet_table(tmp_path, ["PEPTIDEK", "PEPTLDEK"]))) == 2


# --------------------------------------------------------------------------
# End to end: the defect, and the fix, through the CLI
# --------------------------------------------------------------------------
def il_split_inputs(tmp_path: Path) -> list[str]:
    """A Comet I-form and a Casanovo L-form of ONE peptide landing on DIFFERENT regions.

    This is the measured failure: 170 such pairs on the mouse benchmark, of which only
    3 of 340 resulting rows carried two-engine support.
    """
    fasta = tmp_path / "library.fasta"
    fasta.write_text(
        ">sp|P1|A prot one\nAAAPEPTIDEKAAA\n"
        ">sp|P2|B prot two\nCCCPEPTLDEKCCC\n"
        f">{DECOY}sp|P1|A prot one\nAAAKEDITPEPAAA\n"
    )
    diamond = diamond_table(
        tmp_path,
        [("PEPTIDEK", "sp|P1|A", 4, 11, 1e-09), ("PEPTLDEK", "sp|P2|B", 4, 11, 1e-07)],
    )
    return [
        str(comet_table(tmp_path, ["PEPTIDEK"])),
        str(casanovo_table(tmp_path, ["PEPTLDEK"])),
        str(diamond),
        str(fasta),
        DECOY,
    ]


def run_cli(argv: list[str], capsys: pytest.CaptureFixture[str]) -> list[dict[str, str]]:
    assert build_reset_input.main(argv) == 0
    lines = [line for line in capsys.readouterr().out.splitlines() if line]
    return [dict(zip(OUTPUT_COLUMNS, line.split("\t"), strict=True)) for line in lines[1:]]


def test_one_peptide_split_across_two_regions_becomes_one_row(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rows = run_cli(il_split_inputs(tmp_path), capsys)
    assert len(rows) == 1
    row = rows[0]
    assert row["num_comet_peptides"] == "1"
    assert row["num_casanovo_peptides"] == "1", "two-engine support, which the split lost"
    assert row["Peptide"] == "PEPTIDEK", "the better-scoring alignment supplies the region"


def test_the_opt_out_reproduces_the_split(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    rows = run_cli([*il_split_inputs(tmp_path), "--no_il_equivalence"], capsys)
    assert len(rows) == 2
    assert all(
        row["num_comet_peptides"] == "0" or row["num_casanovo_peptides"] == "0" for row in rows
    ), "neither row sees both engines -- this is the defect"


def test_a_pair_landing_on_one_region_was_already_correct(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The 89% case: both spellings already aggregated onto one row, so nothing moves."""
    fasta = tmp_path / "library.fasta"
    fasta.write_text(">sp|P1|A prot one\nAAAPEPTIDEKAAA\n")
    argv = [
        str(comet_table(tmp_path, ["PEPTIDEK"])),
        str(casanovo_table(tmp_path, ["PEPTLDEK"])),
        str(
            diamond_table(
                tmp_path,
                [("PEPTIDEK", "sp|P1|A", 4, 11, 1e-09), ("PEPTLDEK", "sp|P1|A", 4, 11, 1e-07)],
            )
        ),
        str(fasta),
        DECOY,
    ]
    with_flag = run_cli(argv, capsys)
    without = run_cli([*argv, "--no_il_equivalence"], capsys)
    assert len(with_flag) == len(without) == 1
    assert with_flag[0]["num_comet_peptides"] == without[0]["num_comet_peptides"] == "1"
    assert with_flag[0]["num_casanovo_peptides"] == without[0]["num_casanovo_peptides"] == "1"


def test_the_writer_reports_no_missing_peptides_when_spellings_merge(tmp_path: Path) -> None:
    """A regression guard for the join itself: if the tables were normalised but the
    DIAMOND map was not, every I-containing peptide would be reported missing."""
    err = io.StringIO()
    argv = il_split_inputs(tmp_path)
    diamond_map = read_hits(argv[2], DECOY, il_equivalence=True)
    add_subject_sequences(diamond_map, argv[3])
    write_reset_input(
        read_comet_peptides(argv[0], il_equivalence=True),
        read_casanovo_peptides(argv[1], il_equivalence=True),
        diamond_map,
        DECOY,
        out=io.StringIO(),
        err=err,
    )
    assert "missing" not in err.getvalue()
