"""Collect distinct peptide sequences from result tables into a query FASTA."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path


def read_peptide_sequences(file_paths: Iterable[str | Path]) -> set[str]:
    """Union of the first column of each file, header line excluded."""
    sequences: set[str] = set()

    for file_path in file_paths:
        with Path(file_path).open() as handle:
            next(handle)  # header
            for line in handle:
                columns = line.strip().split("\t")
                if columns:
                    sequences.add(columns[0])

    return sequences


def format_fasta(sequences: Iterable[str]) -> Iterator[str]:
    """Emit each peptide as its own entry, named after itself.

    DIAMOND carries the query name through to its output, so naming each entry
    after its sequence is what lets the results be joined back to the peptide.
    """
    for sequence in sequences:
        yield f">{sequence}"
        yield sequence
