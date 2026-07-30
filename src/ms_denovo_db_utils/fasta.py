"""FASTA reading and reversed-decoy generation, transparently gzip-aware."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO, TextIO

GZIP_MAGIC = b"\x1f\x8b"


def is_gzipped(path: str | Path) -> bool:
    """Detect gzip by magic number rather than by file extension."""
    try:
        with Path(path).open("rb") as handle:
            return handle.read(2) == GZIP_MAGIC
    except OSError:
        return False


@contextmanager
def open_text(path: str | Path) -> Iterator[TextIO]:
    """Open a FASTA file for reading, decompressing when needed."""
    handle: TextIO
    # SIM115: this *is* the context manager; the handle is closed below.
    handle = gzip.open(path, "rt") if is_gzipped(path) else Path(path).open()  # noqa: SIM115
    try:
        yield handle
    finally:
        handle.close()


def iter_entries(handle: TextIO) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` pairs; the header keeps its leading '>'.

    Entries with an empty sequence are skipped, matching long-standing
    behaviour of this tool.
    """
    header = ""
    sequence = ""

    for line in handle:
        if line.startswith(">"):
            if header and sequence:
                yield header, sequence
            header = line.strip()
            sequence = ""
        else:
            sequence += line.strip()

    if header and sequence:
        yield header, sequence


def protein_name(header: str) -> str:
    """First whitespace-delimited token of a header, without the '>'."""
    return header[1:].split(maxsplit=1)[0]


def reverse_sequence(sequence: str) -> str:
    return sequence[::-1]


#: Residues trypsin cleaves after. Proline is *not* excepted, matching the
#: Trypsin/P setting the pipeline's comet.params uses.
CLEAVAGE_RESIDUES = "KR"

#: How decoy sequences are derived from their target.
PROTEIN_REVERSAL = "protein"
PEPTIDE_REVERSAL = "peptide"
REVERSAL_METHODS = (PEPTIDE_REVERSAL, PROTEIN_REVERSAL)


def iter_peptides(sequence: str, cleavage_residues: str = CLEAVAGE_RESIDUES) -> Iterator[str]:
    """Split a protein after each cleavage residue, keeping the residue attached."""
    start = 0
    for i, residue in enumerate(sequence):
        if residue in cleavage_residues:
            yield sequence[start : i + 1]
            start = i + 1
    if start < len(sequence):
        yield sequence[start:]


def pseudo_reverse_sequence(sequence: str, cleavage_residues: str = CLEAVAGE_RESIDUES) -> str:
    """Reverse each tryptic peptide in place, holding its last residue fixed.

    This is what makes a decoy peptide the partner of a specific target peptide
    rather than an unrelated string: it preserves length, composition and the
    C-terminal K/R, and — because it works peptide by peptide — a peptide shared
    between two databases yields the same decoy in both, regardless of the
    protein context it appears in. Whole-protein reversal does neither: it moves
    the cleavage boundaries, so the decoys of a shared peptide differ between
    databases and no longer pair one-to-one with targets.

    It also matches how the decoys were built for the training set of the custom
    Casanovo model (`fix_termini="cterm"`), which is trained to be unable to tell
    a peptide from its reversal. That indifference only transfers to the decoys
    the pipeline actually searches if they are made the same way.
    """
    return "".join(
        peptide[:-1][::-1] + peptide[-1] if len(peptide) > 1 else peptide
        for peptide in iter_peptides(sequence, cleavage_residues)
    )


def make_decoy(sequence: str, method: str = PEPTIDE_REVERSAL) -> str:
    if method == PEPTIDE_REVERSAL:
        return pseudo_reverse_sequence(sequence)
    if method == PROTEIN_REVERSAL:
        return reverse_sequence(sequence)
    raise ValueError(f"Unknown reversal method {method!r}; choose from {REVERSAL_METHODS}")


def write_with_decoys(
    input_file: str | Path,
    decoy_prefix: str,
    output: IO[str],
    method: str = PEPTIDE_REVERSAL,
) -> None:
    """Write every target entry followed by its decoy."""
    with open_text(input_file) as handle:
        for header, sequence in iter_entries(handle):
            # A trailing stop character is not a residue and must not become
            # part of the decoy.
            trimmed = sequence.rstrip("*")
            output.write(f"{header}\n{trimmed}\n")
            output.write(f">{decoy_prefix}{header[1:]}\n{make_decoy(trimmed, method)}\n")
