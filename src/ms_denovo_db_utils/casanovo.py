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
    "modifications",
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

#: Casanovo 5.x writes bare residues in ``sequence`` and reports every modification
#: positionally in its own ``modifications`` column:
#: ``1-Oxidation (M):UNIMOD:35; 8-Oxidation (M):UNIMOD:35``.
#:
#: Casanovo 4.x wrote the modification inline in the sequence instead
#: (``M+15.995DLGEEHFK``), and this module used to recover the bare peptide by
#: stripping ``[^A-Z]``. Against 5.x output that strip is a **no-op** -- which is how
#: ``num_peptidoforms`` came to count charge states and nothing else: every peptidoform
#: of one peptide arrived with an identical ``sequence`` field, so the set they were
#: added to could not tell them apart. Measured on the mouse benchmark, 363
#: ``(sequence, charge)`` pairs on the custom checkpoint and 772 on the stock one carried
#: two or more distinct modification states and were counted as one.
#:
#: 5.x is the contract now, so a sequence carrying anything but residues is **rejected
#: rather than stripped**. Stripping is what hid the defect: a 4.x file read by this code
#: would parse without complaint and produce peptidoform counts wrong in exactly the way
#: this replaces.
_NON_RESIDUE = re.compile(r"[^A-Z]")

#: mzTab spells an absent value ``null``.
_MZTAB_NULL = "null"


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


def check_sequence(sequence: str, source: str) -> None:
    """Reject a sequence that is not bare residues.

    See ``_NON_RESIDUE``: an inline modification mass means pre-5.x output, which
    this module no longer reads.
    """
    if _NON_RESIDUE.search(sequence):
        raise MzTabFormatError(
            f"Sequence {sequence!r} in {source} contains a non-residue character. "
            f"This reads Casanovo 5.x mzTab, which writes bare residues in 'sequence' "
            f"and modifications in the 'modifications' column; an inline modification "
            f"mass means pre-5.x output, which is not supported."
        )


def normalise_modifications(modifications: str) -> str:
    """Canonical form of an mzTab ``modifications`` field.

    Two PSMs are the same peptidoform when they carry the same modifications at the
    same positions. Casanovo writes them in position order, so equal peptidoforms
    already produce equal strings; sorting the parts makes the key independent of that
    ordering. Sorting cannot merge two genuinely different peptidoforms, because it
    preserves the multiset of parts.

    ``null`` and an empty field both mean unmodified.
    """
    if modifications.strip() in ("", _MZTAB_NULL):
        return ""
    return ";".join(sorted(part.strip() for part in modifications.split(";") if part.strip()))


def peptidoform_key(sequence: str, modifications: str, charge: int) -> str:
    """Identity of one peptidoform: residues, modification state and charge.

    ``|`` separates the parts because it occurs in none of them. mzTab modification
    syntax uses ``-``, ``;`` and ``:``, so the old ``f"{peptidoform}-{charge}"`` would
    have been ambiguous once the modification column was read.
    """
    return f"{sequence}|{normalise_modifications(modifications)}|{charge}"


def adjust_score(score: float) -> float:
    """Lift negative Casanovo scores into [0, 1].

    Casanovo before 5.2.0 subtracts 1 from the score of a PSM whose peptide
    mass disagrees with the precursor, making it negative. Adding 1 back
    restores the magnitude.

    It does **not** recover what an unpenalised run would have produced, and
    earlier comments here claiming it "effectively disables that filter" were
    wrong. In those versions the precursor check also ends beams early and
    orders the candidate heap; since unpenalised scores lie in [0, 1] and
    penalised ones in [-1, 0], any precursor-matching candidate outranks every
    non-matching one whatever its raw quality. Casanovo has therefore already
    chosen which peptide to report, using penalised scores over a pruned
    search space, before this code sees anything.

    Casanovo 5.2.0 removed the precursor mass filter from de novo mode
    entirely (PR #575), so its de novo scores are never negative and this is a
    no-op on the output this module now reads. Kept as a guard rather than
    deleted alongside the rest of the 4.x handling: it costs one comparison,
    and a negative score reaching the ranking unadjusted would silently invert
    that peptide's rank.
    """
    return score + 1 if score < 0 else score


def process_files(file_paths: Iterable[str | Path]) -> dict[str, CasanovoPeptide]:
    """Read Casanovo mzTab files and return the best PSM per peptide.

    Sorted for the same reason as :func:`~.comet.process_files`: the best-PSM test is a
    strict ``>``, so a tie keeps whichever PSM arrived first, and between files that was
    the caller's argument order. No cross-file tie exists on the mouse benchmark -- both
    orders reproduce the published ``casanovo_peptides.txt`` exactly -- but the exposure
    is real rather than theoretical, because a Casanovo score is a product of per-residue
    probabilities and underflows to exactly 0.0 on a long peptide: 888 of the benchmark's
    113,839. Two spectra in different files predicting one such peptide would tie at 0.0
    and their ``mz_ppm_error`` would be chosen by staging order.
    """
    peptides: dict[str, CasanovoPeptide] = {}

    for file_path in sorted(file_paths, key=str):
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

                sequence = row[index["sequence"]]
                check_sequence(sequence, source)
                modifications = row[index["modifications"]]
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
                current.peptidoforms.add(peptidoform_key(sequence, modifications, charge))

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
