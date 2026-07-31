"""Build the Percolator-style feature table consumed by percolator_RESET.

One output row is one *library* peptide region -- the subject subsequence a
DIAMOND alignment landed on -- rather than one PSM. Comet and Casanovo evidence
for every query peptide that aligned to that region is aggregated onto it, and
the target/decoy label comes from whether the annotated-library protein carries
the decoy prefix. FDR is therefore estimated over annotated-database
identifications, not over spectrum identifications.

Note that the label is a property of the *library* side alone. A query
peptide's own decoy status never affects it, and does not gate whether its
evidence is aggregated: target and decoy rows are assembled by identical rules,
which is what lets decoy rows stand in for null target rows.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from .diamond import DiamondHit

OUTPUT_COLUMNS = (
    "SpecId",
    "Label",
    "ScanNr",
    "database_peptide_length",
    # Both come from the *same* alignment -- the group's highest-bit-score hit --
    # so neither is a maximum taken independently. `max_diamond_perc_identity`
    # was actively untrue: on a group whose alignments run 80% and 99% identity
    # it reported 80%, the identity of the top-scoring hit rather than the top
    # identity. Renamed rather than changed, because features describing one
    # real alignment beat a chimera assembled from several.
    "best_diamond_bitscore",
    "best_diamond_perc_identity",
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
    peptides: Iterable[str],
    diamond_map: Mapping[str, DiamondHit],
) -> dict[str, list[str]]:
    """Map each library subsequence to the query peptides that aligned to it.

    Keys and members are both ordered, so downstream aggregation cannot depend
    on set-iteration order.
    """
    groups: dict[str, list[str]] = {}

    for peptide in sorted(peptides):
        hit = diamond_map.get(peptide)
        if hit is None:
            continue
        if hit.ssequence is None:
            raise ValueError(f"Peptide {peptide} does not have an 'ssequence' property")
        groups.setdefault(hit.ssequence, []).append(peptide)

    return dict(sorted(groups.items()))


def best_diamond_hit(group: Iterable[str], diamond_map: Mapping[str, DiamondHit]) -> DiamondHit:
    """Highest-bit-score alignment in a group, ties broken by subject name."""
    return max(
        (diamond_map[peptide] for peptide in group),
        key=lambda hit: (hit.bitscore, hit.sseqid),
    )


def write_reset_input(  # noqa: PLR0913, PLR0915 - one linear row-assembly routine
    comet_map: Mapping[str, CometRecord],
    casanovo_map: Mapping[str, CasanovoRecord],
    diamond_map: Mapping[str, DiamondHit],
    library_decoy_prefix: str,
    *,
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
        print(
            f"Warning: The following {len(missing_peptides)} peptides are missing "
            f"from diamond_map and will be skipped: " + ", ".join(sorted(missing_peptides)),
            file=err,
        )

    groups = group_by_library_peptide(peptides, diamond_map)

    print("\t".join(OUTPUT_COLUMNS), file=out)
    scan_nr = 1

    # Cross-class evidence: a Comet peptide whose decoy status disagrees with
    # the label of the library region it landed on. Both directions are counted
    # into the features (see the note in the loop); these tallies exist so the
    # frequency is observable, since it is the input to deciding whether the
    # symmetric treatment is the right one.
    decoy_comet_on_target_row = 0
    target_comet_on_decoy_row = 0

    for library_hit_peptide, group in groups.items():
        # Resolved over the whole group before any Comet hit is examined. The
        # decoy filter below reads this, so computing it in the same pass made
        # the result depend on the order the group happened to be iterated in.
        best_hit = best_diamond_hit(group, diamond_map)
        best_diamond_peptide_length = best_hit.subject_length
        best_diamond_bit_score = best_hit.bitscore
        best_diamond_perc_identity = best_hit.pident
        best_diamond_label = -1 if best_hit.is_decoy(library_decoy_prefix) else 1

        # Every distinct protein that produced this region, not just the best
        # hit's. The row's identity is the library peptide, and the same
        # subsequence can occur in several proteins -- common in a
        # metaproteomics library, where related organisms share conserved
        # regions. Naming only one discards annotations the pipeline exists to
        # produce.
        #
        # Taken from the group's *retained* hits rather than from every
        # alignment each query had. That is what makes the list correct by
        # construction: a group member's ssequence is this region by
        # definition, so each protein named genuinely contains it. A query's
        # other alignments may have landed on a different region of a different
        # protein, and listing those would assert something false.
        proteins = sorted({diamond_map[peptide].sseqid for peptide in group})
        best_diamond_proteins = ",".join(proteins)

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

        for peptide in group:
            casanovo_data = casanovo_map.get(peptide)
            comet_data = comet_map.get(peptide)

            if casanovo_data is not None:
                casanovo_num_spectra += casanovo_data.num_spectra
                casanovo_num_peptidoforms += casanovo_data.num_peptidoforms
                casanovo_num_peptides += 1

                if casanovo_data.score > casanovo_best_score:
                    casanovo_best_score = casanovo_data.score
                    casanovo_ppm_error = casanovo_data.mz_ppm_error
                    casanovo_best_rank_score = casanovo_data.rank_score

            if comet_data is not None:
                # Comet evidence is counted whatever its decoy status, matching
                # the transitivity argument in the manuscript: if a spectrum
                # matches a peptide well and that peptide resembles a library
                # peptide, the match is attributable to the latter.
                #
                # This used to drop the decoy-Comet-on-target-row case alone.
                # That was asymmetric -- the mirror case, a target Comet peptide
                # on a library decoy row, was always counted -- so library target
                # and library decoy rows were built from different rules. Decoy
                # rows are the surrogate for null target rows, and that
                # substitution assumes both are assembled the same way.
                if comet_data.is_decoy and best_diamond_label == 1:
                    decoy_comet_on_target_row += 1
                elif not comet_data.is_decoy and best_diamond_label == -1:
                    target_comet_on_decoy_row += 1

                comet_num_spectra += comet_data.num_spectra
                comet_num_peptidoforms += comet_data.num_peptidoforms
                comet_num_peptides += 1

                if comet_data.score > comet_best_score:
                    comet_best_score = comet_data.score
                    comet_ppm_error = comet_data.mz_ppm_error
                    comet_best_rank_score = comet_data.rank_score
                    comet_n_tryptic = comet_data.tryptic_n
                    comet_c_tryptic = comet_data.tryptic_c

        if casanovo_num_spectra == 0 and comet_num_spectra == 0:
            continue

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
            best_diamond_proteins,
        ]

        print("\t".join(str(value) for value in row), file=out)
        scan_nr += 1

    if decoy_comet_on_target_row or target_comet_on_decoy_row:
        print(
            "Note: counted cross-class Comet evidence -- "
            f"{decoy_comet_on_target_row} decoy peptide(s) on library target rows, "
            f"{target_comet_on_decoy_row} target peptide(s) on library decoy rows.",
            file=err,
        )
