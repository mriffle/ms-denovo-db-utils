"""Mass arithmetic and rank scoring shared by the Comet and Casanovo processors."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence

MASS_OF_PROTON = 1.00727647
DELTA_MASS_13C = 1.003355

#: Number of isotope peaks considered when matching a precursor. The
#: monoisotopic peak is often not the most intense one, so the instrument may
#: have selected the first, second or third 13C peak instead.
ISOTOPE_PEAKS_CONSIDERED = 4


def calculate_mz(neutral_mass: float, charge: int) -> float:
    """Convert a neutral mass to m/z at the given charge."""
    return (neutral_mass + (charge * MASS_OF_PROTON)) / charge


def calculate_error_ppm(expected_mz: float, observed_mz: float, charge: int) -> float:
    """Signed precursor mass error in ppm, tolerant of 13C peak misassignment.

    Returns the smallest-magnitude error across the monoisotopic peak and the
    next few isotope peaks, so a precursor selected one 13C too high reports a
    small error rather than a ~1 Da one. The sign of the chosen error is kept.
    """
    min_error_ppm = float("inf")

    for num_13c in range(ISOTOPE_PEAKS_CONSIDERED):
        adjusted_expected_mz = expected_mz + (num_13c * DELTA_MASS_13C) / charge
        error_ppm = (observed_mz - adjusted_expected_mz) / adjusted_expected_mz * 1000000

        if abs(error_ppm) < abs(min_error_ppm):
            min_error_ppm = error_ppm

    return min_error_ppm


def rank_scores(values: Sequence[float], *, higher_is_better: bool) -> dict[float, float]:
    """Map each distinct score to its competition rank, normalised by count.

    Ties share the best rank and consume the ranks below them (1, 2, 2, 4), and
    the result is divided by the number of values so it lands in (0, 1].
    """
    if not values:
        return {}

    occurrences = Counter(values)
    total = len(values)

    ranks: dict[float, float] = {}
    current_rank = 1
    for value in sorted(occurrences, reverse=higher_is_better):
        ranks[value] = current_rank / total
        current_rank += occurrences[value]

    return ranks
