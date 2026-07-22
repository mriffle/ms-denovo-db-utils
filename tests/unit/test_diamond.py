"""DIAMOND result parsing, best-hit selection and ambiguity filtering."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from ms_denovo_db_utils.diamond import (
    DiamondHit,
    add_subject_sequences,
    parse_line,
    read_hits,
)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
GOLDEN = FIXTURES / "golden"
DECOY = "LIBRARY_DECOY_"


def hit_line(
    query: str,
    subject: str,
    *,
    evalue: float,
    bitscore: float,
    sstart: int = 1,
    send: int = 14,
    pident: float = 100.0,
) -> str:
    return "\t".join(
        [
            query,
            subject,
            str(pident),
            "14",
            "0",
            "0",
            "1",
            "14",
            str(sstart),
            str(send),
            repr(evalue),
            str(bitscore),
        ]
    )


def write_hits(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "hits.dmnd.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------
def test_parse_line_reads_the_columns_diamond_actually_emits() -> None:
    hit = parse_line(hit_line("PEPTIDEK", "sp|P1|X", evalue=1e-9, bitscore=60.0), 1)
    assert (hit.qseqid, hit.sseqid) == ("PEPTIDEK", "sp|P1|X")
    assert (hit.sstart, hit.send, hit.evalue, hit.bitscore) == (1, 14, 1e-9, 60.0)


def test_parse_line_rejects_the_wrong_column_count() -> None:
    with pytest.raises(ValueError, match="Expected 12 columns, but found 3"):
        parse_line("a\tb\tc", 7)


def test_subject_length_is_inclusive_of_both_endpoints() -> None:
    hit = parse_line(hit_line("P", "sp|X", evalue=1.0, bitscore=1.0, sstart=10, send=20), 1)
    assert hit.subject_length == 11


def test_comments_and_blank_lines_are_skipped(tmp_path: Path) -> None:
    path = write_hits(
        tmp_path,
        ["# a comment", "", hit_line("PEPTIDEK", "sp|P1|X", evalue=1e-9, bitscore=60.0)],
    )
    assert set(read_hits(path, DECOY)) == {"PEPTIDEK"}


# --------------------------------------------------------------------------
# Best-hit selection
# --------------------------------------------------------------------------
def test_lowest_evalue_wins(tmp_path: Path) -> None:
    path = write_hits(
        tmp_path,
        [
            hit_line("PEPTIDEK", "sp|WORSE|X", evalue=1e-3, bitscore=90.0),
            hit_line("PEPTIDEK", "sp|BETTER|X", evalue=1e-9, bitscore=30.0),
        ],
    )
    assert read_hits(path, DECOY)["PEPTIDEK"].sseqid == "sp|BETTER|X"


@pytest.mark.parametrize("reverse", [False, True])
def test_tied_evalues_resolve_the_same_way_in_either_input_order(
    tmp_path: Path, reverse: bool
) -> None:
    lines = [
        hit_line("PEPTIDEK", "sp|A|X", evalue=1e-9, bitscore=50.0),
        hit_line("PEPTIDEK", "sp|B|X", evalue=1e-9, bitscore=70.0),
    ]
    path = write_hits(tmp_path, list(reversed(lines)) if reverse else lines)
    # Same e-value, so the higher bit score decides -- not the input order.
    assert read_hits(path, DECOY)["PEPTIDEK"].sseqid == "sp|B|X"


@pytest.mark.parametrize("reverse", [False, True])
def test_fully_tied_hits_resolve_the_same_way_in_either_input_order(
    tmp_path: Path, reverse: bool
) -> None:
    lines = [
        hit_line("PEPTIDEK", "sp|B|X", evalue=1e-9, bitscore=50.0),
        hit_line("PEPTIDEK", "sp|A|X", evalue=1e-9, bitscore=50.0),
    ]
    path = write_hits(tmp_path, list(reversed(lines)) if reverse else lines)
    assert read_hits(path, DECOY)["PEPTIDEK"].sseqid == "sp|A|X"


# --------------------------------------------------------------------------
# Ambiguity filtering -- the defect fixed alongside these tests
# --------------------------------------------------------------------------
@pytest.mark.parametrize("reverse", [False, True])
def test_peptide_hitting_both_target_and_decoy_is_dropped_in_either_order(
    tmp_path: Path, reverse: bool
) -> None:
    """A peptide matching both carries no target/decoy signal, either way round.

    Previously the accumulated protein set was cleared whenever a better hit
    arrived, so this was only caught when the rows happened to arrive with the
    better hit first.
    """
    lines = [
        hit_line("PEPTIDEK", "sp|P1|X", evalue=1e-3, bitscore=40.0),
        hit_line("PEPTIDEK", f"{DECOY}sp|P1|X", evalue=1e-9, bitscore=60.0),
    ]
    path = write_hits(tmp_path, list(reversed(lines)) if reverse else lines)
    assert "PEPTIDEK" not in read_hits(path, DECOY)


def test_both_orders_of_the_same_ambiguity_are_treated_alike() -> None:
    """The regression fixture: two peptides, mirror-image row orders."""
    kept = read_hits(FIXTURES / "diamond/ambiguous.dmnd.txt", DECOY)
    assert kept == {}


def test_peptide_hitting_only_decoys_is_kept(tmp_path: Path) -> None:
    path = write_hits(
        tmp_path,
        [
            hit_line("PEPTIDEK", f"{DECOY}sp|A|X", evalue=1e-9, bitscore=60.0),
            hit_line("PEPTIDEK", f"{DECOY}sp|B|X", evalue=1e-5, bitscore=40.0),
        ],
    )
    assert "PEPTIDEK" in read_hits(path, DECOY)


def test_peptide_hitting_only_targets_is_kept(tmp_path: Path) -> None:
    path = write_hits(
        tmp_path,
        [
            hit_line("PEPTIDEK", "sp|A|X", evalue=1e-9, bitscore=60.0),
            hit_line("PEPTIDEK", "sp|B|X", evalue=1e-5, bitscore=40.0),
        ],
    )
    assert "PEPTIDEK" in read_hits(path, DECOY)


def test_ambiguity_of_one_peptide_does_not_affect_another(tmp_path: Path) -> None:
    path = write_hits(
        tmp_path,
        [
            hit_line("AMBIGK", "sp|A|X", evalue=1e-3, bitscore=40.0),
            hit_line("AMBIGK", f"{DECOY}sp|A|X", evalue=1e-9, bitscore=60.0),
            hit_line("CLEANK", "sp|B|X", evalue=1e-9, bitscore=60.0),
        ],
    )
    assert set(read_hits(path, DECOY)) == {"CLEANK"}


# --------------------------------------------------------------------------
# Subject sequence resolution
# --------------------------------------------------------------------------
def test_subject_sequence_is_the_aligned_region_of_the_protein() -> None:
    hits = read_hits(FIXTURES / "diamond/hits.dmnd.txt", DECOY)
    add_subject_sequences(hits, GOLDEN / "library_plusdecoys.fasta")
    # SAMPLEPEPTIDER sits at residues 12-25 of sp|P1|PROT1.
    assert hits["SAMPLEPEPTIDER"].ssequence == "SAMPLEPEPTIDER"


def test_distinct_queries_can_resolve_to_the_same_subject_region() -> None:
    """This is what makes two query peptides share one output row."""
    hits = read_hits(FIXTURES / "diamond/hits.dmnd.txt", DECOY)
    add_subject_sequences(hits, GOLDEN / "library_plusdecoys.fasta")
    assert hits["SHAREDPEPTIDEK"].ssequence == hits["TIEDPEPTIDEAK"].ssequence


def test_subject_sequences_resolve_from_a_gzipped_library(tmp_path: Path) -> None:
    gz = tmp_path / "library.fasta.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write((GOLDEN / "library_plusdecoys.fasta").read_text())

    hits = read_hits(FIXTURES / "diamond/hits.dmnd.txt", DECOY)
    add_subject_sequences(hits, gz)
    assert hits["SAMPLEPEPTIDER"].ssequence == "SAMPLEPEPTIDER"


def test_hit_on_a_protein_absent_from_the_fasta_is_an_error(tmp_path: Path) -> None:
    """Silently emitting a row with no subject sequence would corrupt output."""
    fasta = tmp_path / "library.fasta"
    fasta.write_text(">sp|OTHER|X desc\nMKVLAAGIVGL\n")
    hits = {"PEPTIDEK": DiamondHit("PEPTIDEK", "sp|MISSING|X", 100.0, 1, 8, 1e-9, 60.0)}
    with pytest.raises(ValueError, match="does not have an 'ssequence'"):
        add_subject_sequences(hits, fasta)
