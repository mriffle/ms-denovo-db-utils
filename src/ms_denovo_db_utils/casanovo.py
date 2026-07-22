"""Collapse Casanovo mzTab output to one best row per distinct peptide."""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from .massutil import calculate_error_ppm, rank_scores

REQUIRED_COLUMNS = (
    "sequence",
    "charge",
    "search_engine_score[1]",
    "calc_mass_to_charge",
    "exp_mass_to_charge",
)

OUTPUT_COLUMNS = (
    "peptide_sequence",
    "charge",
    "search_engine_score[1]",
    "file",
    "mz_ppm_error",
    "num_spectra",
    "rank_score",
    "num_peptidoforms",
)

#: Modifications appear inline in the mzTab sequence, e.g. "M+15.995DLGEEHFK".
_NON_RESIDUE = re.compile(r"[^A-Z]")


class MzTabFormatError(ValueError):
    """Raised when an mzTab file is missing its PSM header or a needed column."""


@dataclass
class CasanovoPeptide:
    """Best PSM for one peptide, plus aggregates over all of its PSMs."""

    charge: int
    score: float
    source_file: str
    mz_ppm_error: float
    num_spectra: int = 0
    peptidoforms: set[str] = field(default_factory=set)
    rank_score: float = 0.0


def plain_sequence(peptidoform: str) -> str:
    """Strip inline modification masses, leaving bare residues."""
    return _NON_RESIDUE.sub("", peptidoform)


def adjust_score(score: float) -> float:
    """Lift negative Casanovo scores into [0, 1].

    Casanovo subtracts 1 from the score when a PSM fails its precursor m/z
    filter. Adding 1 back effectively disables that filter, which is what the
    pipeline wants: the homology search decides what is credible, not Casanovo.
    """
    return score + 1 if score < 0 else score


def process_files(file_paths: Iterable[str | Path]) -> dict[str, CasanovoPeptide]:
    """Read Casanovo mzTab files and return the best PSM per peptide."""
    peptides: dict[str, CasanovoPeptide] = {}

    for file_path in file_paths:
        source = str(file_path)
        with Path(file_path).open() as handle:
            reader = csv.reader(handle, delimiter="\t")

            headers: list[str] | None = None
            for row in reader:
                if row and row[0].startswith("PSH"):
                    headers = row
                    break
            if headers is None:
                raise MzTabFormatError(f"No PSH header line found in {source}")

            try:
                index = {name: headers.index(name) for name in REQUIRED_COLUMNS}
            except ValueError as exc:
                raise MzTabFormatError(f"Missing expected column in {source}: {exc}") from exc

            for row in reader:
                if not row or not row[0].startswith("PSM"):
                    continue

                peptidoform = row[index["sequence"]]
                sequence = plain_sequence(peptidoform)
                charge = int(float(row[index["charge"]]))
                # The ppm error is computed from the reported score's PSM before
                # any adjustment; adjust_score only affects ranking.
                mz_ppm_error = calculate_error_ppm(
                    float(row[index["calc_mass_to_charge"]]),
                    float(row[index["exp_mass_to_charge"]]),
                    charge,
                )
                score = adjust_score(float(row[index["search_engine_score[1]"]]))

                existing = peptides.get(sequence)
                if existing is None or score > existing.score:
                    best = CasanovoPeptide(
                        charge=charge,
                        score=score,
                        source_file=source,
                        mz_ppm_error=mz_ppm_error,
                    )
                    if existing is not None:
                        best.num_spectra = existing.num_spectra
                        best.peptidoforms = existing.peptidoforms
                    peptides[sequence] = best

                current = peptides[sequence]
                current.num_spectra += 1
                current.peptidoforms.add(f"{peptidoform}-{charge}")

    _assign_rank_scores(peptides)
    return peptides


def _assign_rank_scores(peptides: dict[str, CasanovoPeptide]) -> None:
    """Rank peptides by score, best (highest) first."""
    ranks = rank_scores([p.score for p in peptides.values()], higher_is_better=True)
    for peptide in peptides.values():
        peptide.rank_score = ranks[peptide.score]


def format_results(peptides: dict[str, CasanovoPeptide]) -> Iterator[str]:
    """Yield the tab-delimited output lines, header first."""
    yield "\t".join(OUTPUT_COLUMNS)
    for sequence, data in peptides.items():
        yield "\t".join(
            [
                sequence,
                str(data.charge),
                str(data.score),
                data.source_file,
                f"{data.mz_ppm_error:.2f}",
                str(data.num_spectra),
                str(data.rank_score),
                str(len(data.peptidoforms)),
            ]
        )
