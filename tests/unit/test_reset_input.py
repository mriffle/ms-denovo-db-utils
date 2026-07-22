"""Assembly of the percolator_RESET feature table.

The unit of a row is a library peptide region, so most of what matters here is
how evidence from several query peptides is aggregated onto one row and how the
target/decoy label is decided.
"""

from __future__ import annotations

import io
import math
from pathlib import Path

import pytest

from ms_denovo_db_utils.diamond import DiamondHit
from ms_denovo_db_utils.reset_input import (
    OUTPUT_COLUMNS,
    CasanovoRecord,
    CometRecord,
    best_diamond_hit,
    group_by_library_peptide,
    read_casanovo_peptides,
    read_comet_peptides,
    write_reset_input,
)

DECOY = "LIBRARY_DECOY_"

COLUMN = {name: i for i, name in enumerate(OUTPUT_COLUMNS)}


def hit(
    query: str,
    subject: str,
    *,
    ssequence: str,
    bitscore: float = 50.0,
    evalue: float = 1e-9,
    pident: float = 100.0,
    sstart: int = 1,
) -> DiamondHit:
    return DiamondHit(
        qseqid=query,
        sseqid=subject,
        pident=pident,
        sstart=sstart,
        send=sstart + len(ssequence) - 1,
        evalue=evalue,
        bitscore=bitscore,
        ssequence=ssequence,
    )


def comet(
    *,
    e_value: float = 1e-9,
    num_spectra: int = 1,
    num_peptidoforms: int = 1,
    is_decoy: bool = False,
    rank_score: float = 0.5,
) -> CometRecord:
    return CometRecord(
        e_value=e_value,
        num_spectra=num_spectra,
        num_peptidoforms=num_peptidoforms,
        is_decoy=is_decoy,
        rank_score=rank_score,
        mz_ppm_error="1.00",
        tryptic_n="1",
        tryptic_c="1",
    )


def casanovo(
    *,
    score: float = 0.9,
    num_spectra: int = 1,
    num_peptidoforms: int = 1,
    rank_score: float = 0.5,
) -> CasanovoRecord:
    return CasanovoRecord(
        score=score,
        num_spectra=num_spectra,
        num_peptidoforms=num_peptidoforms,
        rank_score=rank_score,
        mz_ppm_error="2.00",
    )


def build(
    comet_map: dict[str, CometRecord],
    casanovo_map: dict[str, CasanovoRecord],
    diamond_map: dict[str, DiamondHit],
) -> tuple[list[dict[str, str]], str]:
    """Run the writer and return parsed rows plus stderr."""
    out, err = io.StringIO(), io.StringIO()
    write_reset_input(comet_map, casanovo_map, diamond_map, DECOY, out=out, err=err)
    lines = [line for line in out.getvalue().splitlines() if line]
    assert lines[0].split("\t") == list(OUTPUT_COLUMNS)
    rows = [dict(zip(OUTPUT_COLUMNS, line.split("\t"), strict=True)) for line in lines[1:]]
    return rows, err.getvalue()


# --------------------------------------------------------------------------
# Comet score transform
# --------------------------------------------------------------------------
def test_smaller_evalues_transform_to_larger_scores() -> None:
    assert comet(e_value=1e-9).score > comet(e_value=1e-5).score > comet(e_value=1e-1).score


def test_an_evalue_of_zero_is_finite() -> None:
    """Comet reports 0.0; the floor keeps the transform from dividing by zero."""
    score = comet(e_value=0.0).score
    assert math.isfinite(score)
    assert score == pytest.approx(20.0)


# --------------------------------------------------------------------------
# Grouping
# --------------------------------------------------------------------------
def test_queries_sharing_a_subject_region_form_one_group() -> None:
    diamond_map = {
        "AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK"),
        "BBBK": hit("BBBK", "sp|P1", ssequence="LIBPEPK"),
        "CCCK": hit("CCCK", "sp|P2", ssequence="OTHERK"),
    }
    groups = group_by_library_peptide(["AAAK", "BBBK", "CCCK"], diamond_map)
    assert groups == {"LIBPEPK": ["AAAK", "BBBK"], "OTHERK": ["CCCK"]}


def test_peptides_without_a_diamond_hit_are_skipped() -> None:
    diamond_map = {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK")}
    assert group_by_library_peptide(["AAAK", "MISSINGK"], diamond_map) == {"LIBPEPK": ["AAAK"]}


def test_an_unresolved_subject_sequence_is_an_error() -> None:
    unresolved = DiamondHit("AAAK", "sp|P1", 100.0, 1, 7, 1e-9, 50.0, ssequence=None)
    with pytest.raises(ValueError, match="does not have an 'ssequence'"):
        group_by_library_peptide(["AAAK"], {"AAAK": unresolved})


def test_best_hit_in_a_group_is_the_highest_bit_score() -> None:
    diamond_map = {
        "AAAK": hit("AAAK", "sp|LOW", ssequence="LIBPEPK", bitscore=30.0),
        "BBBK": hit("BBBK", "sp|HIGH", ssequence="LIBPEPK", bitscore=70.0),
    }
    assert best_diamond_hit(["AAAK", "BBBK"], diamond_map).sseqid == "sp|HIGH"


# --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------
def test_a_target_library_hit_labels_the_row_positive() -> None:
    rows, _ = build({"AAAK": comet()}, {}, {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK")})
    assert rows[0]["Label"] == "1"


def test_a_decoy_library_hit_labels_the_row_negative() -> None:
    rows, _ = build(
        {"AAAK": comet()}, {}, {"AAAK": hit("AAAK", f"{DECOY}sp|P1", ssequence="LIBPEPK")}
    )
    assert rows[0]["Label"] == "-1"


def test_the_label_follows_the_best_hit_not_the_first() -> None:
    """A weaker decoy alignment must not flip a target row's label."""
    diamond_map = {
        "AAAK": hit("AAAK", f"{DECOY}sp|P1", ssequence="LIBPEPK", bitscore=30.0),
        "BBBK": hit("BBBK", "sp|P2", ssequence="LIBPEPK", bitscore=70.0),
    }
    rows, _ = build({"AAAK": comet(), "BBBK": comet()}, {}, diamond_map)
    assert rows[0]["Label"] == "1"
    assert rows[0]["Proteins"] == "sp|P2"


# --------------------------------------------------------------------------
# Aggregation
# --------------------------------------------------------------------------
def test_evidence_from_every_query_in_a_group_is_summed() -> None:
    diamond_map = {
        "AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK"),
        "BBBK": hit("BBBK", "sp|P1", ssequence="LIBPEPK"),
    }
    comet_map = {
        "AAAK": comet(num_spectra=3, num_peptidoforms=2),
        "BBBK": comet(num_spectra=4, num_peptidoforms=1),
    }
    rows, _ = build(comet_map, {}, diamond_map)
    assert len(rows) == 1
    assert rows[0]["comet_num_spectra"] == "7"
    assert rows[0]["comet_num_peptidoforms"] == "3"
    assert rows[0]["num_comet_peptides"] == "2"


def test_casanovo_and_comet_evidence_land_on_the_same_row() -> None:
    diamond_map = {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK")}
    rows, _ = build({"AAAK": comet(num_spectra=2)}, {"AAAK": casanovo(num_spectra=5)}, diamond_map)
    assert rows[0]["comet_num_spectra"] == "2"
    assert rows[0]["casanovo_num_spectra"] == "5"


def test_a_casanovo_only_peptide_produces_a_row_with_no_comet_evidence() -> None:
    diamond_map = {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK")}
    rows, _ = build({}, {"AAAK": casanovo(num_spectra=3)}, diamond_map)
    assert rows[0]["comet_num_spectra"] == "0"
    assert rows[0]["num_comet_peptides"] == "0"
    assert rows[0]["casanovo_num_spectra"] == "3"


def test_the_best_scoring_hit_supplies_the_reported_features() -> None:
    diamond_map = {
        "AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK"),
        "BBBK": hit("BBBK", "sp|P1", ssequence="LIBPEPK"),
    }
    weak = CometRecord(1e-2, 1, 1, False, 0.9, "9.99", "0", "0")
    strong = CometRecord(1e-12, 1, 1, False, 0.1, "0.11", "1", "1")
    rows, _ = build({"AAAK": weak, "BBBK": strong}, {}, diamond_map)
    assert rows[0]["comet_ppm_error"] == "0.11"
    assert rows[0]["comet_n_tryptic"] == "1"


def test_database_peptide_length_is_the_aligned_region_length() -> None:
    diamond_map = {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPTIDEK", sstart=5)}
    rows, _ = build({"AAAK": comet()}, {}, diamond_map)
    assert rows[0]["database_peptide_length"] == "11"


# --------------------------------------------------------------------------
# Decoy Comet hits
# --------------------------------------------------------------------------
def test_a_decoy_comet_hit_on_a_target_row_is_discarded() -> None:
    diamond_map = {
        "AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK", bitscore=70.0),
        "DDDK": hit("DDDK", "sp|P1", ssequence="LIBPEPK", bitscore=30.0),
    }
    comet_map = {"AAAK": comet(num_spectra=2), "DDDK": comet(num_spectra=9, is_decoy=True)}
    rows, err = build(comet_map, {}, diamond_map)
    assert rows[0]["comet_num_spectra"] == "2"
    assert rows[0]["num_comet_peptides"] == "1"
    assert "Ignoring decoy comet hit for DDDK" in err


def test_a_decoy_comet_hit_on_a_decoy_row_is_kept() -> None:
    """Decoy-on-decoy is the consistent case and carries real signal."""
    diamond_map = {"DDDK": hit("DDDK", f"{DECOY}sp|P1", ssequence="LIBPEPK")}
    rows, err = build({"DDDK": comet(num_spectra=9, is_decoy=True)}, {}, diamond_map)
    assert rows[0]["comet_num_spectra"] == "9"
    assert "Ignoring" not in err


def test_a_row_left_with_no_evidence_is_dropped_entirely() -> None:
    """Discarding the only Comet hit must not leave a zero-evidence row."""
    diamond_map = {"DDDK": hit("DDDK", "sp|P1", ssequence="LIBPEPK")}
    rows, err = build({"DDDK": comet(is_decoy=True)}, {}, diamond_map)
    assert rows == []
    assert "Ignoring decoy comet hit" in err


# --------------------------------------------------------------------------
# Combined rank score and row ordering
# --------------------------------------------------------------------------
def test_combined_rank_score_rewards_good_ranks_from_both_engines() -> None:
    diamond_map = {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK")}
    rows, _ = build(
        {"AAAK": comet(rank_score=0.1)}, {"AAAK": casanovo(rank_score=0.2)}, diamond_map
    )
    assert float(rows[0]["combined_rank_score"]) == pytest.approx(4 - 0.1 - 0.2)


def test_a_missing_engine_contributes_the_worst_rank() -> None:
    diamond_map = {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK")}
    rows, _ = build({"AAAK": comet(rank_score=0.25)}, {}, diamond_map)
    assert float(rows[0]["combined_rank_score"]) == pytest.approx(4 - 0.25 - 2)


def test_rows_are_ordered_by_spec_id_with_matching_scan_numbers() -> None:
    diamond_map = {
        "AAAK": hit("AAAK", "sp|P1", ssequence="ZZZPEPK"),
        "BBBK": hit("BBBK", "sp|P2", ssequence="AAAPEPK"),
        "CCCK": hit("CCCK", "sp|P3", ssequence="MMMPEPK"),
    }
    comet_map = {name: comet() for name in diamond_map}
    rows, _ = build(comet_map, {}, diamond_map)
    assert [row["SpecId"] for row in rows] == ["AAAPEPK_1", "MMMPEPK_1", "ZZZPEPK_1"]
    assert [row["ScanNr"] for row in rows] == ["1", "2", "3"]


def test_missing_peptides_are_reported_sorted() -> None:
    diamond_map = {"AAAK": hit("AAAK", "sp|P1", ssequence="LIBPEPK")}
    comet_map = {"AAAK": comet(), "ZZZK": comet(), "MMMK": comet()}
    _, err = build(comet_map, {}, diamond_map)
    assert "2 peptides are missing" in err
    assert "MMMK, ZZZK" in err


# --------------------------------------------------------------------------
# Reading the intermediate tables
# --------------------------------------------------------------------------
def test_comet_table_round_trips_through_the_reader(tmp_path: Path) -> None:
    path = tmp_path / "comet_peptides.txt"
    path.write_text(
        "plain_peptide\tcharge\te-value\tprotein\tfile\ttryptic_n\ttryptic_c\tnum_spectra\t"
        "mz_ppm_error\tis_decoy\tproteins\trank_score\tnum_peptidoforms\n"
        "PEPTIDEK\t2\t1e-09\tsp|A\trun.txt\t1\t0\t4\t-1.50\t1\tsp|A\t0.25\t3\n"
    )
    record = read_comet_peptides(path)["PEPTIDEK"]
    assert record.e_value == 1e-9
    assert record.num_spectra == 4
    assert record.num_peptidoforms == 3
    assert record.is_decoy is True
    assert record.rank_score == 0.25
    assert record.mz_ppm_error == "-1.50"
    assert (record.tryptic_n, record.tryptic_c) == ("1", "0")


def test_casanovo_table_round_trips_through_the_reader(tmp_path: Path) -> None:
    path = tmp_path / "casanovo_peptides.txt"
    path.write_text(
        "peptide_sequence\tcharge\tsearch_engine_score[1]\tfile\tmz_ppm_error\t"
        "num_spectra\trank_score\tnum_peptidoforms\n"
        "PEPTIDEK\t2\t0.9905\trun.mztab\t0.90\t4\t0.4\t3\n"
    )
    record = read_casanovo_peptides(path)["PEPTIDEK"]
    assert record.score == 0.9905
    assert record.num_spectra == 4
    assert record.rank_score == 0.4
    assert record.mz_ppm_error == "0.90"
