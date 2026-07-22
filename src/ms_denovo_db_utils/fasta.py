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


def write_with_decoys(
    input_file: str | Path,
    decoy_prefix: str,
    output: IO[str],
) -> None:
    """Write every target entry followed by its reversed decoy."""
    with open_text(input_file) as handle:
        for header, sequence in iter_entries(handle):
            # A trailing stop character is not a residue and must not become
            # the first residue of the decoy.
            trimmed = sequence.rstrip("*")
            output.write(f"{header}\n{trimmed}\n")
            output.write(f">{decoy_prefix}{header[1:]}\n{reverse_sequence(trimmed)}\n")
