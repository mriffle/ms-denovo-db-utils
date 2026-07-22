"""FASTA parsing and reversed-decoy generation, including gzip handling."""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import pytest

from ms_denovo_db_utils.fasta import (
    is_gzipped,
    iter_entries,
    open_text,
    protein_name,
    reverse_sequence,
    write_with_decoys,
)

DECOY_PREFIX = "LIBRARY_DECOY_"


def decoys_for(path: Path, prefix: str = DECOY_PREFIX) -> str:
    buffer = io.StringIO()
    write_with_decoys(path, prefix, buffer)
    return buffer.getvalue()


# --------------------------------------------------------------------------
# Gzip detection
# --------------------------------------------------------------------------
def test_gzip_is_detected_by_magic_number_not_extension(tmp_path: Path) -> None:
    """The pipeline passes files whose names do not always reflect encoding."""
    misnamed = tmp_path / "not_obviously_gzipped.fasta"
    with gzip.open(misnamed, "wt") as handle:
        handle.write(">sp|A|X\nPEPTIDEK\n")
    assert is_gzipped(misnamed)


def test_a_plain_file_is_not_reported_as_gzipped(tmp_path: Path) -> None:
    path = tmp_path / "plain.fasta"
    path.write_text(">sp|A|X\nPEPTIDEK\n")
    assert not is_gzipped(path)


def test_a_missing_file_is_not_reported_as_gzipped(tmp_path: Path) -> None:
    assert not is_gzipped(tmp_path / "absent.fasta")


def test_gzipped_and_plain_files_read_identically(tmp_path: Path) -> None:
    text = ">sp|A|X desc\nPEPTIDEK\n>sp|B|Y\nMKVLAAGIVGL\n"
    plain = tmp_path / "plain.fasta"
    plain.write_text(text)
    gzipped = tmp_path / "gzipped.fasta.gz"
    with gzip.open(gzipped, "wt") as handle:
        handle.write(text)

    with open_text(plain) as a, open_text(gzipped) as b:
        assert list(iter_entries(a)) == list(iter_entries(b))


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------
def test_sequences_spanning_multiple_lines_are_joined() -> None:
    handle = io.StringIO(">sp|A|X\nPEPT\nIDEK\nMKVL\n")
    assert list(iter_entries(handle)) == [(">sp|A|X", "PEPTIDEKMKVL")]


def test_the_final_entry_is_emitted() -> None:
    handle = io.StringIO(">sp|A|X\nPEPTIDEK\n>sp|B|Y\nMKVLAAGIVGL\n")
    assert len(list(iter_entries(handle))) == 2


def test_an_entry_with_no_sequence_is_skipped() -> None:
    handle = io.StringIO(">sp|EMPTY|X\n>sp|B|Y\nPEPTIDEK\n")
    assert list(iter_entries(handle)) == [(">sp|B|Y", "PEPTIDEK")]


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        (">sp|P1|PROT1 Test protein one", "sp|P1|PROT1"),
        (">sp|P1|PROT1", "sp|P1|PROT1"),
        (">tr|Q9X|X_HUMAN Some long description here", "tr|Q9X|X_HUMAN"),
    ],
)
def test_protein_name_is_the_first_token(header: str, expected: str) -> None:
    """DIAMOND reports this same token as sseqid, which is what joins them."""
    assert protein_name(header) == expected


def test_reverse_sequence_reverses() -> None:
    assert reverse_sequence("PEPTIDEK") == "KEDITPEP"


# --------------------------------------------------------------------------
# Decoy generation
# --------------------------------------------------------------------------
def test_each_target_is_followed_by_its_decoy(tmp_path: Path) -> None:
    path = tmp_path / "in.fasta"
    path.write_text(">sp|A|X desc\nPEPTIDEK\n")
    assert decoys_for(path) == (f">sp|A|X desc\nPEPTIDEK\n>{DECOY_PREFIX}sp|A|X desc\nKEDITPEP\n")


def test_the_decoy_keeps_the_full_description(tmp_path: Path) -> None:
    """DIAMOND takes sseqid from the first token, which must stay unique."""
    path = tmp_path / "in.fasta"
    path.write_text(">sp|A|X a long description\nPEPTIDEK\n")
    assert f">{DECOY_PREFIX}sp|A|X a long description" in decoys_for(path)


def test_a_trailing_stop_codon_is_stripped_before_reversal(tmp_path: Path) -> None:
    """Otherwise '*' becomes the decoy's first residue and DIAMOND sees junk."""
    path = tmp_path / "in.fasta"
    path.write_text(">sp|A|X\nPEPTIDEK*\n")
    output = decoys_for(path)
    assert "PEPTIDEK\n" in output
    assert "KEDITPEP\n" in output
    assert "*" not in output


def test_every_protein_gets_a_decoy(tmp_path: Path) -> None:
    path = tmp_path / "in.fasta"
    path.write_text(">sp|A|X\nPEPTIDEK\n>sp|B|Y\nMKVLAAGIVGL\n>sp|C|Z\nAAGGK\n")
    output = decoys_for(path)
    assert output.count(">") == 6
    assert output.count(f">{DECOY_PREFIX}") == 3


def test_the_decoy_prefix_is_configurable(tmp_path: Path) -> None:
    path = tmp_path / "in.fasta"
    path.write_text(">sp|A|X\nPEPTIDEK\n")
    assert ">rev_sp|A|X" in decoys_for(path, "rev_")


def test_a_gzipped_input_is_read_transparently(tmp_path: Path) -> None:
    path = tmp_path / "in.fasta.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">sp|A|X\nPEPTIDEK\n")
    assert "KEDITPEP" in decoys_for(path)


def test_reversal_is_an_involution_on_the_decoy(tmp_path: Path) -> None:
    """A decoy of a decoy is the original, so the reversal is not lossy."""
    sequence = "MKVLAAGIVGLSAMPLEPEPTIDER"
    assert reverse_sequence(reverse_sequence(sequence)) == sequence
