"""DIAMOND result parsing, best-hit selection and ambiguity filtering."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from ms_denovo_db_utils.diamond import (
    DECOY as CLASS_DECOY,
)
from ms_denovo_db_utils.diamond import (
    ENTRAPMENT_DECOY as CLASS_ENTRAPMENT_DECOY,
)
from ms_denovo_db_utils.diamond import (
    ENTRAPMENT_TARGET as CLASS_ENTRAPMENT_TARGET,
)
from ms_denovo_db_utils.diamond import (
    MIXED as CLASS_MIXED,
)
from ms_denovo_db_utils.diamond import (
    TARGET as CLASS_TARGET,
)
from ms_denovo_db_utils.diamond import (
    DiamondHit,
    add_subject_sequences,
    classify_subject,
    parse_line,
    read_hits,
    region_entrapment_axis,
    region_entrapment_class,
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


# --------------------------------------------------------------------------
# Entrapment classification (ENTRAPMENT_DESIGN.md §2.3, §2.5)
# --------------------------------------------------------------------------
LIB_DECOY = DECOY  # the prefix; CLASS_DECOY above is the class name
ENT = "ENTRAPMENT_"


def test_classify_subject_covers_the_four_classes() -> None:
    """The prefixes compose with the decoy prefix OUTERMOST, because entrapments are
    injected before decoy generation and the decoy step wraps them."""
    assert classify_subject("sp|P1|F", LIB_DECOY, ENT) == CLASS_TARGET
    assert classify_subject(f"{ENT}sp|P1|F", LIB_DECOY, ENT) == CLASS_ENTRAPMENT_TARGET
    assert classify_subject(f"{LIB_DECOY}sp|P1|F", LIB_DECOY, ENT) == CLASS_DECOY
    assert classify_subject(f"{LIB_DECOY}{ENT}sp|P1|F", LIB_DECOY, ENT) == CLASS_ENTRAPMENT_DECOY


def test_classify_subject_strips_the_decoy_prefix_before_testing_entrapment() -> None:
    """The mistake the design warns about. Testing ENTRAPMENT_ against a raw accession
    puts every decoyed entrapment in the plain-decoy bucket, so N_E loses its own decoys
    and the 1:1 correspondence RESET's estimate rests on is silently broken."""
    decoyed = f"{LIB_DECOY}{ENT}sp|P1|F"
    assert not decoyed.startswith(ENT)
    assert classify_subject(decoyed, LIB_DECOY, ENT) == CLASS_ENTRAPMENT_DECOY


def test_classify_subject_without_an_entrapment_prefix_yields_only_two_classes() -> None:
    assert classify_subject(f"{ENT}sp|P1|F", LIB_DECOY) == CLASS_TARGET
    assert classify_subject(f"{LIB_DECOY}{ENT}sp|P1|F", LIB_DECOY) == CLASS_DECOY


def test_an_entrapment_target_is_a_target_for_the_label() -> None:
    """Entrapments must be indistinguishable to RESET -- that is the whole point. They
    compete as targets; entrapment membership is a separate axis, recovered downstream."""
    for name in (f"{ENT}sp|P1|F", "sp|P1|F"):
        assert classify_subject(name, LIB_DECOY, ENT) in (CLASS_TARGET, CLASS_ENTRAPMENT_TARGET)


def test_a_region_from_one_class_keeps_that_class() -> None:
    got = region_entrapment_class([f"{ENT}sp|P1|F", f"{ENT}sp|P2|G"], LIB_DECOY, ENT)
    assert got == CLASS_ENTRAPMENT_TARGET


def test_a_region_reachable_from_both_an_original_and_an_entrapment_is_mixed() -> None:
    """§2.5. Such a region is not cleanly attributable and must not be counted as an
    entrapment discovery -- nor as a clean target one."""
    got = region_entrapment_class(["sp|P1|F", f"{ENT}sp|P2|G"], LIB_DECOY, ENT)
    assert got == CLASS_MIXED


def test_a_region_naming_no_protein_raises() -> None:
    with pytest.raises(ValueError, match="at least one protein"):
        region_entrapment_class([], LIB_DECOY, ENT)


def test_a_target_decoy_region_is_not_reported_as_entrapment_mixed() -> None:
    """The bug this caught. A region CAN name both a target and a decoy protein -- two
    query peptides, each individually unambiguous, landing on the same subsequence in
    both. It is pre-existing and documented (SPECIFICATION.md §5), and measured at 2 of
    196,085 regions on a library with no entrapments in it at all. Counting it as `mixed`
    would pollute the one number that decides whether the experiment is compromised."""
    proteins = ["sp|P1|F", f"{LIB_DECOY}sp|P2|G"]
    assert region_entrapment_axis(proteins, LIB_DECOY, ENT) == "original"
    assert region_entrapment_class(proteins, LIB_DECOY, ENT) == CLASS_TARGET
    assert region_entrapment_class(proteins, LIB_DECOY, ENT, is_decoy_row=True) == CLASS_DECOY


def test_the_reported_class_takes_its_label_from_the_row_not_the_protein_list() -> None:
    """`Label` is decided by the best hit alone. Re-deriving it from the protein list here
    could disagree with the row the reader is looking at."""
    proteins = [f"{ENT}sp|P1|F"]
    assert region_entrapment_class(proteins, LIB_DECOY, ENT) == CLASS_ENTRAPMENT_TARGET
    assert (
        region_entrapment_class(proteins, LIB_DECOY, ENT, is_decoy_row=True)
        == CLASS_ENTRAPMENT_DECOY
    )


def test_the_entrapment_axis_ignores_decoy_status_entirely() -> None:
    mixed_labels_one_axis = [f"{ENT}sp|P1|F", f"{LIB_DECOY}{ENT}sp|P2|G"]
    assert region_entrapment_axis(mixed_labels_one_axis, LIB_DECOY, ENT) == "entrapment"


def test_the_drop_rule_is_unchanged_by_the_entrapment_prefix(tmp_path: Path) -> None:
    """The property step 3 of the design turns on: default-off must be bit-identical, and
    passing the prefix must not change WHICH peptides are retained. Entrapment only ever
    splits a class into two that fall on the same side of the target/decoy test."""
    lines = [
        # spans target/decoy -> dropped, with and without the prefix
        ("Q1", "sp|P1|F"),
        ("Q1", f"{LIB_DECOY}sp|P2|G"),
        # spans only the entrapment axis within targets -> retained
        ("Q2", "sp|P3|H"),
        ("Q2", f"{ENT}sp|P4|I"),
        # spans only the entrapment axis within decoys -> retained
        ("Q3", f"{LIB_DECOY}sp|P5|J"),
        ("Q3", f"{LIB_DECOY}{ENT}sp|P6|K"),
        # entrapment target vs real decoy -> still spans the axis, dropped
        ("Q4", f"{ENT}sp|P7|L"),
        ("Q4", f"{LIB_DECOY}sp|P8|M"),
    ]
    path = tmp_path / "hits.dmnd.txt"
    with path.open("w") as handle:
        for qseqid, sseqid in lines:
            handle.write(f"{qseqid}\t{sseqid}\t90.0\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50.0\n")

    without = read_hits(path, LIB_DECOY)
    with_prefix = read_hits(path, LIB_DECOY, entrapment_prefix=ENT)
    assert set(without) == {"Q2", "Q3"}
    assert set(with_prefix) == set(without)
    assert {q: h.sseqid for q, h in with_prefix.items()} == {
        q: h.sseqid for q, h in without.items()
    }
