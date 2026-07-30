"""Collapse Comet PSM output to one best row per distinct peptide.

Comet reports one row per PSM, ranked within each spectrum. Only each
spectrum's best-scoring match is considered (see :data:`RANK_COLUMN`); of
those, this keeps the lowest-e-value PSM for each peptide and attaches counts
and features derived from every PSM that peptide appeared in.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

from .massutil import calculate_error_ppm, calculate_mz, rank_scores

REQUIRED_COLUMNS = (
    "plain_peptide",
    "charge",
    "e-value",
    "protein",
    "modified_peptide",
    "calc_neutral_mass",
    "exp_neutral_mass",
)

#: Comet's within-spectrum rank; rank 1 is that spectrum's best-scoring match.
#: Deliberately *not* in REQUIRED_COLUMNS -- a result file written without it
#: still parses, it just cannot be filtered. Comet decides this ranking (by
#: xcorr), and taking its word for it is what makes "top-scoring PSM per
#: spectrum" mean the same thing here as it does in Comet.
RANK_COLUMN = "num"
TOP_RANK = 1

OUTPUT_COLUMNS = (
    "plain_peptide",
    "charge",
    "e-value",
    "protein",
    "file",
    "tryptic_n",
    "tryptic_c",
    "num_spectra",
    "mz_ppm_error",
    "is_decoy",
    "proteins",
    "rank_score",
    "num_peptidoforms",
)


class CometFormatError(ValueError):
    """Raised when a Comet result file lacks a column this tool needs."""


@dataclass
class CometPeptide:
    """Best PSM for one peptide, plus aggregates over all of its PSMs."""

    charge: int
    e_value: float
    protein: str
    source_file: str
    tryptic_n: int
    tryptic_c: int
    mz_ppm_error: float
    is_decoy: int
    num_spectra: int = 0
    peptidoforms: set[str] = field(default_factory=set)
    rank_score: float = 0.0


def is_n_tryptic(modified_peptide: str) -> int:
    """Whether the residue preceding the peptide implies a tryptic N-terminus.

    Comet's ``modified_peptide`` carries the flanking residues, as in
    ``K.PEPTIDER.A``; ``-`` marks a protein terminus, which counts as tryptic.
    """
    return 1 if modified_peptide.startswith(("R", "K", "-")) else 0


def is_c_tryptic(plain_peptide: str, modified_peptide: str) -> int:
    """Whether the peptide ends in R/K, or at the protein C-terminus."""
    return 1 if (plain_peptide.endswith(("R", "K")) or modified_peptide.endswith("-")) else 0


def is_decoy(protein: str, decoy_prefix: str) -> int:
    """A peptide counts as a decoy only when every protein it maps to is one."""
    proteins = protein.split(",")
    return int(all(p.startswith(decoy_prefix) for p in proteins))


def _column_indices(headers: Sequence[str], source: str) -> dict[str, int]:
    try:
        return {name: headers.index(name) for name in REQUIRED_COLUMNS}
    except ValueError as exc:
        raise CometFormatError(f"Missing expected column in {source}: {exc}") from exc


def is_top_ranked(row: Sequence[str], rank_index: int | None, source: str) -> bool:
    """Whether this row is its spectrum's best-scoring match.

    Rows ranked below the best are further matches to a spectrum already
    accounted for; counting them would inflate ``num_spectra`` (which would no
    longer count spectra) and ``num_peptidoforms``. Ties at rank 1 are all
    kept, since each is genuinely a top-scoring match.

    Returns True when the file carries no rank column, leaving such files to
    behave as they always have.
    """
    if rank_index is None:
        return True
    try:
        return int(row[rank_index]) == TOP_RANK
    except ValueError as exc:
        raise CometFormatError(
            f"Non-numeric {RANK_COLUMN!r} value in {source}: {row[rank_index]!r}"
        ) from exc


def _iter_rows(handle: TextIO) -> Iterator[list[str]]:
    """Yield data rows, skipping Comet's version banner and header line."""
    next(handle)  # "CometVersion ..." banner
    reader = csv.reader(handle, delimiter="\t")
    headers = next(reader)
    yield headers
    yield from reader


def process_files(file_paths: Iterable[str | Path], decoy_prefix: str) -> dict[str, CometPeptide]:
    """Read Comet result files and return the best PSM per peptide."""
    peptides: dict[str, CometPeptide] = {}

    for file_path in file_paths:
        source = str(file_path)
        with Path(file_path).open() as handle:
            rows = _iter_rows(handle)
            headers = next(rows)
            index = _column_indices(headers, source)
            rank_index = headers.index(RANK_COLUMN) if RANK_COLUMN in headers else None

            for row in rows:
                if not is_top_ranked(row, rank_index, source):
                    continue

                plain_peptide = row[index["plain_peptide"]]
                charge = int(row[index["charge"]])
                e_value = float(row[index["e-value"]])
                protein = row[index["protein"]]
                modified_peptide = row[index["modified_peptide"]]
                calc_mz = calculate_mz(float(row[index["calc_neutral_mass"]]), charge)
                exp_mz = calculate_mz(float(row[index["exp_neutral_mass"]]), charge)

                existing = peptides.get(plain_peptide)
                if existing is None or e_value < existing.e_value:
                    best = CometPeptide(
                        charge=charge,
                        e_value=e_value,
                        protein=protein,
                        source_file=source,
                        tryptic_n=is_n_tryptic(modified_peptide),
                        tryptic_c=is_c_tryptic(plain_peptide, modified_peptide),
                        mz_ppm_error=calculate_error_ppm(calc_mz, exp_mz, charge),
                        is_decoy=is_decoy(protein, decoy_prefix),
                    )
                    if existing is not None:
                        best.num_spectra = existing.num_spectra
                        best.peptidoforms = existing.peptidoforms
                    peptides[plain_peptide] = best

                current = peptides[plain_peptide]
                current.num_spectra += 1
                current.peptidoforms.add(f"{modified_peptide}-{charge}")

    _assign_rank_scores(peptides)
    return peptides


def _assign_rank_scores(peptides: dict[str, CometPeptide]) -> None:
    """Rank peptides by e-value, best (lowest) first."""
    ranks = rank_scores([p.e_value for p in peptides.values()], higher_is_better=False)
    for peptide in peptides.values():
        peptide.rank_score = ranks[peptide.e_value]


def format_results(peptides: dict[str, CometPeptide]) -> Iterator[str]:
    """Yield the tab-delimited output lines, header first."""
    yield "\t".join(OUTPUT_COLUMNS)
    for sequence, data in peptides.items():
        yield "\t".join(
            [
                sequence,
                str(data.charge),
                str(data.e_value),
                data.protein,
                data.source_file,
                str(data.tryptic_n),
                str(data.tryptic_c),
                str(data.num_spectra),
                f"{data.mz_ppm_error:.2f}",
                str(data.is_decoy),
                # Deliberately repeats `protein`. Kept because comet_peptides.txt
                # is published to results/ and consumers outside this repo may
                # read it; removing a column is a breaking format change.
                data.protein,
                str(data.rank_score),
                str(len(data.peptidoforms)),
            ]
        )
