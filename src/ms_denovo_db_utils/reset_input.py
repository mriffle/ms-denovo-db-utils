"""Build the Percolator-style feature table consumed by percolator_RESET.

One output row is one *library* peptide region -- the subject subsequence a
DIAMOND alignment landed on -- rather than one PSM. Comet and Casanovo evidence
for every query peptide that aligned to that region is aggregated onto it, and
the target/decoy label comes from whether the annotated-library protein carries
the decoy prefix. FDR is therefore estimated over annotated-database
identifications, not over spectrum identifications.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from .diamond import DiamondHit

OUTPUT_COLUMNS = (
    "SpecId",
    "Label",
    "ScanNr",
    "database_peptide_length",
    "max_diamond_bitscore",
    "max_diamond_perc_identity",
    "num_casanovo_peptides",
    "num_comet_peptides",
    "casanovo_num_spectra",
    "casanovo_best_score",
    "casanovo_ppm_error",
    "casanovo_num_peptidoforms",
    "comet_num_spectra",
    "comet_n_tryptic",
    "comet_c_tryptic",
    "comet_best_score",
    "comet_ppm_error",
    "comet_num_peptidoforms",
    "combined_rank_score",
    "Peptide",
    "Proteins",
)

#: Guards against division by zero for Comet e-values of exactly 0, which the
#: search engine does report.
_EVALUE_FLOOR = 1e-20

#: Worst possible rank_score, used when an engine contributed no hit at all.
_WORST_RANK = 2.0


@dataclass(frozen=True)
class CometRecord:
    """One row of comet_peptides.txt."""

    e_value: float
    num_spectra: int
    num_peptidoforms: int
    is_decoy: bool
    rank_score: float
    # Written straight through to the output, so kept as text to preserve the
    # upstream formatting exactly.
    mz_ppm_error: str
    tryptic_n: str
    tryptic_c: str

    @property
    def score(self) -> float:
        """Comet e-values transformed so that larger is better."""
        return math.log10(1 + (1 / (self.e_value + _EVALUE_FLOOR)))


@dataclass(frozen=True)
class CasanovoRecord:
    """One row of casanovo_peptides.txt."""

    score: float
    num_spectra: int
    num_peptidoforms: int
    rank_score: float
    mz_ppm_error: str


def _read_table(file_path: str | Path) -> Iterator[tuple[str, dict[str, str]]]:
    with Path(file_path).open() as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            columns = line.strip().split("\t")
            yield columns[0], {header[i]: columns[i] for i in range(1, len(header))}


def read_comet_peptides(file_path: str | Path) -> dict[str, CometRecord]:
    return {
        peptide: CometRecord(
            e_value=float(row["e-value"]),
            num_spectra=int(row["num_spectra"]),
            num_peptidoforms=int(row["num_peptidoforms"]),
            is_decoy=row["is_decoy"] == "1",
            rank_score=float(row["rank_score"]),
            mz_ppm_error=row["mz_ppm_error"],
            tryptic_n=row["tryptic_n"],
            tryptic_c=row["tryptic_c"],
        )
        for peptide, row in _read_table(file_path)
    }


def read_casanovo_peptides(file_path: str | Path) -> dict[str, CasanovoRecord]:
    return {
        peptide: CasanovoRecord(
            score=float(row["search_engine_score[1]"]),
            num_spectra=int(row["num_spectra"]),
            num_peptidoforms=int(row["num_peptidoforms"]),
            rank_score=float(row["rank_score"]),
            mz_ppm_error=row["mz_ppm_error"],
        )
        for peptide, row in _read_table(file_path)
    }


def group_by_library_peptide(
    peptides: Iterator[str] | set[str],
    diamond_map: Mapping[str, DiamondHit],
) -> dict[str, set[str]]:
    """Map each library subsequence to the query peptides that aligned to it."""
    groups: dict[str, set[str]] = {}

    for peptide in peptides:
        hit = diamond_map.get(peptide)
        if hit is None:
            continue
        if hit.ssequence is None:
            raise ValueError(f"Peptide {peptide} does not have an 'ssequence' property")
        groups.setdefault(hit.ssequence, set()).add(peptide)

    return groups


def write_reset_input(  # noqa: PLR0913, PLR0915 - one linear row-assembly routine
    comet_map: Mapping[str, CometRecord],
    casanovo_map: Mapping[str, CasanovoRecord],
    diamond_map: Mapping[str, DiamondHit],
    library_decoy_prefix: str,
    out: TextIO | None = None,
    err: TextIO | None = None,
) -> None:
    # Resolved here rather than as default arguments, which would bind
    # whatever sys.stdout happened to be when this module was imported.
    out = sys.stdout if out is None else out
    err = sys.stderr if err is None else err

    peptides = set(comet_map.keys()) | set(casanovo_map.keys())
    missing_peptides = peptides - set(diamond_map.keys())

    if missing_peptides:
        # BUG(#2): set iteration order, so this listing differs between runs.
        print(
            f"Warning: The following {len(missing_peptides)} peptides are missing "
            f"from diamond_map and will be skipped: " + ", ".join(missing_peptides),
            file=err,
        )

    groups = group_by_library_peptide(peptides, diamond_map)

    print("\t".join(OUTPUT_COLUMNS), file=out)
    scan_nr = 1

    # BUG(#2): `groups` is keyed in set-iteration order and each group is itself
    # a set, so row order, ScanNr, and -- because best_diamond_label is both
    # computed and read inside this loop -- the decoy filtering below all vary
    # between runs. Fixed in the next commit.
    for library_hit_peptide, group in groups.items():
        best_diamond_peptide_length = 0
        best_diamond_bit_score: float = 0
        best_diamond_perc_identity: float = 0
        best_diamond_label = 1
        best_diamond_protein = ""

        casanovo_num_spectra = 0
        casanovo_num_peptides = 0
        casanovo_best_score: float = 0
        casanovo_ppm_error = "0"
        casanovo_num_peptidoforms = 0
        casanovo_best_rank_score = _WORST_RANK

        comet_num_spectra = 0
        comet_num_peptides = 0
        comet_n_tryptic = "0"
        comet_c_tryptic = "0"
        comet_best_score: float = 0
        comet_ppm_error = "0"
        comet_num_peptidoforms = 0
        comet_best_rank_score = _WORST_RANK
        comet_best_is_decoy = False

        for peptide in group:
            diamond_hit = diamond_map[peptide]
            casanovo_data = casanovo_map.get(peptide)
            comet_data = comet_map.get(peptide)

            if diamond_hit.bitscore > best_diamond_bit_score:
                best_diamond_bit_score = diamond_hit.bitscore
                best_diamond_peptide_length = diamond_hit.subject_length
                best_diamond_perc_identity = diamond_hit.pident
                best_diamond_label = -1 if diamond_hit.is_decoy(library_decoy_prefix) else 1
                best_diamond_protein = diamond_hit.sseqid

            if casanovo_data is not None:
                casanovo_num_spectra += casanovo_data.num_spectra
                casanovo_num_peptidoforms += casanovo_data.num_peptidoforms
                casanovo_num_peptides += 1

                if casanovo_data.score > casanovo_best_score:
                    casanovo_best_score = casanovo_data.score
                    casanovo_ppm_error = casanovo_data.mz_ppm_error
                    casanovo_best_rank_score = casanovo_data.rank_score

            if comet_data is not None:
                if comet_data.is_decoy and best_diamond_label == 1:
                    print(
                        f"Warning: Ignoring decoy comet hit for {peptide} that "
                        f"matched target in library: {best_diamond_protein}",
                        file=err,
                    )
                else:
                    comet_num_spectra += comet_data.num_spectra
                    comet_num_peptidoforms += comet_data.num_peptidoforms
                    comet_num_peptides += 1

                    if comet_data.score > comet_best_score:
                        comet_best_score = comet_data.score
                        comet_ppm_error = comet_data.mz_ppm_error
                        comet_best_rank_score = comet_data.rank_score
                        comet_n_tryptic = comet_data.tryptic_n
                        comet_c_tryptic = comet_data.tryptic_c
                        comet_best_is_decoy = comet_data.is_decoy

        if casanovo_num_spectra == 0 and comet_num_spectra == 0:
            continue

        if comet_best_is_decoy and best_diamond_label == 1:
            raise ValueError(
                "Found situation where best comet hit is a decoy and best "
                "diamond hit is a target, stopping."
            )

        combined_rank_score = 4 - comet_best_rank_score - casanovo_best_rank_score

        row = [
            f"{library_hit_peptide}_1",
            best_diamond_label,
            scan_nr,
            best_diamond_peptide_length,
            best_diamond_bit_score,
            best_diamond_perc_identity,
            casanovo_num_peptides,
            comet_num_peptides,
            casanovo_num_spectra,
            casanovo_best_score,
            casanovo_ppm_error,
            casanovo_num_peptidoforms,
            comet_num_spectra,
            comet_n_tryptic,
            comet_c_tryptic,
            comet_best_score,
            comet_ppm_error,
            comet_num_peptidoforms,
            combined_rank_score,
            library_hit_peptide,
            best_diamond_protein,
        ]

        print("\t".join(str(value) for value in row), file=out)
        scan_nr += 1
