"""Guards against reintroducing quadratic behaviour on realistic input sizes.

Rank scoring once counted occurrences by scanning the full value list for every
distinct value. On a search producing hundreds of thousands of peptides that is
quadratic; these bound it loosely enough not to flake on a busy CI runner but
tightly enough that a return to scanning would fail outright.
"""

from __future__ import annotations

import time

import pytest

from ms_denovo_db_utils.massutil import rank_scores

pytestmark = pytest.mark.slow

PEPTIDE_COUNT = 200_000
TIME_BUDGET_SECONDS = 10.0


def test_rank_scoring_scales_to_a_realistic_peptide_count() -> None:
    values = [float(i) for i in range(PEPTIDE_COUNT)]

    started = time.perf_counter()
    ranks = rank_scores(values, higher_is_better=False)
    elapsed = time.perf_counter() - started

    assert len(ranks) == PEPTIDE_COUNT
    assert elapsed < TIME_BUDGET_SECONDS, f"rank scoring took {elapsed:.1f}s"


def test_rank_scoring_scales_when_almost_every_value_is_tied() -> None:
    """Heavy ties are the worst case for a per-value rescan."""
    values = [float(i % 10) for i in range(PEPTIDE_COUNT)]

    started = time.perf_counter()
    ranks = rank_scores(values, higher_is_better=True)
    elapsed = time.perf_counter() - started

    assert len(ranks) == 10
    assert elapsed < TIME_BUDGET_SECONDS, f"rank scoring took {elapsed:.1f}s"
