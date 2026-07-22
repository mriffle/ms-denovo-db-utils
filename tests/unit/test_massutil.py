"""Mass arithmetic and rank scoring.

The ppm calculation is the piece with real behaviour in it: it deliberately
tolerates the instrument having selected an isotope peak rather than the
monoisotopic one, and that tolerance has to work without hiding genuine error.
"""

from __future__ import annotations

import math

import pytest

from ms_denovo_db_utils.massutil import (
    DELTA_MASS_13C,
    MASS_OF_PROTON,
    calculate_error_ppm,
    calculate_mz,
    rank_scores,
)


# --------------------------------------------------------------------------
# m/z
# --------------------------------------------------------------------------
@pytest.mark.parametrize("charge", [1, 2, 3, 4])
def test_mz_adds_one_proton_per_charge(charge: int) -> None:
    neutral = 1000.0
    assert calculate_mz(neutral, charge) == pytest.approx(
        (neutral + charge * MASS_OF_PROTON) / charge
    )


def test_singly_charged_mz_is_the_neutral_mass_plus_a_proton() -> None:
    assert calculate_mz(1000.0, 1) == pytest.approx(1000.0 + MASS_OF_PROTON)


# --------------------------------------------------------------------------
# ppm error
# --------------------------------------------------------------------------
def test_exact_match_has_no_error() -> None:
    assert calculate_error_ppm(500.0, 500.0, 2) == pytest.approx(0.0)


@pytest.mark.parametrize("ppm", [-10.0, -5.0, -0.5, 0.5, 5.0, 10.0])
def test_error_is_recovered_with_its_sign(ppm: float) -> None:
    expected_mz = 500.0
    observed_mz = expected_mz * (1 + ppm / 1e6)
    assert calculate_error_ppm(expected_mz, observed_mz, 2) == pytest.approx(ppm, abs=1e-6)


@pytest.mark.parametrize("num_13c", [1, 2, 3])
@pytest.mark.parametrize("charge", [1, 2, 3])
def test_isotope_peak_selection_is_forgiven(num_13c: int, charge: int) -> None:
    """A precursor one to three 13C high must report a near-zero error."""
    expected_mz = 500.0
    observed_mz = expected_mz + num_13c * DELTA_MASS_13C / charge
    assert calculate_error_ppm(expected_mz, observed_mz, charge) == pytest.approx(0.0, abs=1e-6)


def test_a_fourth_isotope_peak_is_not_forgiven() -> None:
    """Only the monoisotopic peak plus three isotopes are considered.

    Beyond that the mismatch is real and must show up as a large error rather
    than being silently absorbed.
    """
    charge = 2
    observed_mz = 500.0 + 4 * DELTA_MASS_13C / charge
    assert abs(calculate_error_ppm(500.0, observed_mz, charge)) > 100


def test_isotope_tolerance_does_not_mask_a_genuine_offset() -> None:
    """A peptide 0.5 Da off matches no isotope peak and stays a large error."""
    assert abs(calculate_error_ppm(500.0, 500.25, 2)) > 100


def test_smallest_magnitude_error_wins_across_isotopes() -> None:
    """Slightly below one 13C should report a small negative error."""
    charge = 2
    expected_mz = 500.0
    observed_mz = expected_mz + DELTA_MASS_13C / charge - 500.0 * 3e-6
    result = calculate_error_ppm(expected_mz, observed_mz, charge)
    assert result == pytest.approx(-3.0, abs=0.1)


# --------------------------------------------------------------------------
# rank scoring
# --------------------------------------------------------------------------
def test_empty_input_produces_no_ranks() -> None:
    assert rank_scores([], higher_is_better=True) == {}


def test_lower_is_better_ranks_ascending() -> None:
    ranks = rank_scores([1e-9, 1e-5, 1e-1], higher_is_better=False)
    assert ranks[1e-9] < ranks[1e-5] < ranks[1e-1]


def test_higher_is_better_ranks_descending() -> None:
    ranks = rank_scores([0.9, 0.5, 0.1], higher_is_better=True)
    assert ranks[0.9] < ranks[0.5] < ranks[0.1]


def test_ties_share_a_rank_and_consume_the_ranks_below_them() -> None:
    """Competition ranking: 1, 2, 2, 4 -- not 1, 2, 2, 3."""
    ranks = rank_scores([10.0, 20.0, 20.0, 40.0], higher_is_better=False)
    assert ranks[10.0] == pytest.approx(1 / 4)
    assert ranks[20.0] == pytest.approx(2 / 4)
    assert ranks[40.0] == pytest.approx(4 / 4)


def test_scores_are_normalised_into_the_unit_interval() -> None:
    values = [float(i) for i in range(10)]
    ranks = rank_scores(values, higher_is_better=True)
    assert all(0 < rank <= 1 for rank in ranks.values())
    assert max(ranks.values()) == pytest.approx(1.0)


def test_a_single_value_ranks_at_one() -> None:
    assert rank_scores([0.5], higher_is_better=True) == {0.5: 1.0}


def test_all_values_tied_share_the_best_rank() -> None:
    ranks = rank_scores([7.0, 7.0, 7.0, 7.0], higher_is_better=False)
    assert ranks == {7.0: pytest.approx(0.25)}


def test_rank_is_position_over_count_not_a_percentile() -> None:
    ranks = rank_scores([1.0, 2.0, 3.0, 4.0, 5.0], higher_is_better=False)
    assert [ranks[v] for v in (1.0, 2.0, 3.0, 4.0, 5.0)] == [
        pytest.approx(x) for x in (0.2, 0.4, 0.6, 0.8, 1.0)
    ]


def test_infinity_and_zero_are_ordered_sensibly() -> None:
    ranks = rank_scores([0.0, math.inf], higher_is_better=False)
    assert ranks[0.0] < ranks[math.inf]
