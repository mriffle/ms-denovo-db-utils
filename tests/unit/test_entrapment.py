"""Entrapment generation: shuffling, determinism, uniqueness and the pairing record.

The properties tested here are the ones the FDR-validation experiment rests on. An
entrapment that is not verifiably false, or that differs between two databases, does not
fail loudly -- it produces a plausible-looking FDP that is wrong.
"""

from __future__ import annotations

import gzip
import io
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from ms_denovo_db_utils.cli import generate_entrapments
from ms_denovo_db_utils.entrapment import (
    PAIR_TAG,
    EntrapmentStats,
    _record_verbatim_runs,
    collect_peptides,
    entrapment_header,
    make_entrapment,
    shuffle_peptide,
    write_with_entrapments,
)
from ms_denovo_db_utils.fasta import (
    iter_entries,
    iter_peptides,
    make_decoy,
    open_text,
    protein_name,
)

ENTRAPMENT_PREFIX = "ENTRAPMENT_"


def entrapments_for(path: Path, **kwargs: object) -> str:
    buffer = io.StringIO()
    kwargs.setdefault("err", None)
    write_with_entrapments(path, ENTRAPMENT_PREFIX, buffer, **kwargs)  # type: ignore[arg-type]
    return buffer.getvalue()


def write_fasta(path: Path, entries: dict[str, str]) -> Path:
    with path.open("w") as handle:
        for name, sequence in entries.items():
            handle.write(f">{name}\n{sequence}\n")
    return path


# --------------------------------------------------------------------------
# shuffle_peptide
# --------------------------------------------------------------------------
def test_shuffle_holds_the_c_terminal_residue_fixed() -> None:
    """The tryptic K/R must survive, or the entrapment digests differently from the
    original and the two are no longer comparable."""
    for peptide in ("PEPTIDEK", "SAMPLERR", "ACDEFGHIK"):
        assert shuffle_peptide(peptide)[-1] == peptide[-1]


def test_shuffle_preserves_composition() -> None:
    """A permutation, not a substitution: same residues, same length, same mass."""
    peptide = "ACDEFGHIKLMNPQR"
    shuffled = shuffle_peptide(peptide)
    assert Counter(shuffled) == Counter(peptide)
    assert len(shuffled) == len(peptide)


def test_shuffle_actually_moves_something() -> None:
    peptide = "ACDEFGHIKLMNPQR"
    assert shuffle_peptide(peptide) != peptide


def test_shuffle_is_deterministic_for_a_given_seed_and_attempt() -> None:
    assert shuffle_peptide("ACDEFGHIK", 7, 0) == shuffle_peptide("ACDEFGHIK", 7, 0)


def test_shuffle_differs_between_seeds_and_between_attempts() -> None:
    """Retries must explore, and two databases built with different seeds must differ."""
    assert shuffle_peptide("ACDEFGHIKLMNPQR", 0, 0) != shuffle_peptide("ACDEFGHIKLMNPQR", 1, 0)
    assert shuffle_peptide("ACDEFGHIKLMNPQR", 0, 0) != shuffle_peptide("ACDEFGHIKLMNPQR", 0, 1)


@pytest.mark.parametrize("peptide", ["", "K", "AK"])
def test_shuffle_returns_short_peptides_unchanged(peptide: str) -> None:
    """With the C-terminus pinned, a two-residue peptide has one movable position and
    nowhere to move it. Returned unchanged rather than pretending to shuffle."""
    assert shuffle_peptide(peptide) == peptide


def test_shuffle_does_not_depend_on_the_interpreters_hash_seed() -> None:
    """hash() is salted per process; hashlib is not. If this module ever switched to
    hash(), every run would build a different database from the same inputs -- and the
    Comet-side and library-side entrapments for a shared peptide would disagree."""
    code = (
        "from ms_denovo_db_utils.entrapment import shuffle_peptide;"
        "print(shuffle_peptide('ACDEFGHIKLMNPQR', 3, 0))"
    )
    outputs = set()
    for seed in ("0", "1", "12345"):
        proc = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, "PYTHONHASHSEED": seed},
        )
        outputs.add(proc.stdout.strip())
    assert len(outputs) == 1, f"shuffle varies with PYTHONHASHSEED: {outputs}"


# --------------------------------------------------------------------------
# make_entrapment
# --------------------------------------------------------------------------
def test_entrapment_preserves_every_cleavage_site() -> None:
    """Peptide-by-peptide shuffling must not move a K or R, so the entrapment digests
    into the same number of peptides, of the same lengths, as the original."""
    sequence = "ACDEFGHIKLMNPQRSTVWYK" * 3
    got = make_entrapment(sequence)
    assert [len(p) for p in iter_peptides(got)] == [len(p) for p in iter_peptides(sequence)]
    assert Counter(got) == Counter(sequence)


def test_entrapment_is_not_the_reversed_decoy() -> None:
    """The one way to make this experiment vacuous: generate entrapments the same way
    decoys are generated, so every entrapment is byte-identical to a decoy."""
    sequence = "ACDEFGHIKLMNPQRSTVWYK" * 3
    assert make_entrapment(sequence) != make_decoy(sequence)


def test_entrapment_of_the_same_sequence_is_identical_across_calls() -> None:
    sequence = "ACDEFGHIKLMNPQRSTVWYK" * 2
    assert make_entrapment(sequence, seed=5) == make_entrapment(sequence, seed=5)


def test_a_shared_peptide_gets_the_same_entrapment_in_a_different_protein() -> None:
    """The cross-database property. A peptide shared between the annotated library and
    the Comet database must yield the same entrapment in both, or the two databases
    disagree about the same peptide. Protein context must not enter the shuffle."""
    shared = "ACDEFGHILMNPQTVWYK"
    first = make_entrapment("MSTARTK" + shared + "WYVTENDR")
    second = make_entrapment("QQQQQQR" + shared)
    common = set(iter_peptides(first)) & set(iter_peptides(second))
    assert common, "the shared peptide produced no common entrapment peptide"


def test_entrapment_shuffles_peptides_outside_the_window_by_default() -> None:
    """The deviation from the paper, and the reason for it: this pipeline reaches the
    library by homology, so a verbatim stretch is an alignment anchor that lets a real
    query peptide align to an entrapment protein."""
    long_peptide = "ACDEFGHILMNPQTVWYACDEFGHILMNPQTVWYACDEFGHILMNPQTVWY"  # 51 aa, no K/R
    assert len(long_peptide) > 35
    assert make_entrapment(long_peptide) != long_peptide


def test_entrapment_keeps_out_of_window_peptides_when_asked() -> None:
    """The paper's literal construction stays reachable behind a flag."""
    long_peptide = "ACDEFGHILMNPQTVWYACDEFGHILMNPQTVWYACDEFGHILMNPQTVWY"
    got = make_entrapment(long_peptide, shuffle_outside_window=False)
    assert got == long_peptide


def test_keeping_out_of_window_peptides_leaves_far_more_intact(tmp_path: Path) -> None:
    """Pins the measurement that drove the default. Under the paper's window a long
    tryptic peptide survives verbatim; shuffling everything removes it."""
    sequence = "MK" + "ACDEFGHILMNPQTVWY" * 4 + "K"
    kept = EntrapmentStats()
    make_entrapment(sequence, stats=kept, shuffle_outside_window=False)
    shuffled = EntrapmentStats()
    make_entrapment(sequence, stats=shuffled, shuffle_outside_window=True)
    assert kept.residues_intact > shuffled.residues_intact


# --------------------------------------------------------------------------
# Uniqueness
# --------------------------------------------------------------------------
def test_an_entrapment_peptide_never_equals_a_real_peptide(tmp_path: Path) -> None:
    """The property N_E depends on. If an entrapment peptide coincides with a real one,
    a true identification is counted as a false entrapment discovery."""
    entries = {f"sp|P{i:05d}|TEST": "ACDEFGHIKLMNPQRSTVWYK" for i in range(20)}
    entries["sp|P99999|MIX"] = "STVWYACDEFGHIKLMNPQRK"
    path = write_fasta(tmp_path / "in.fasta", entries)
    originals = collect_peptides(path)

    text = entrapments_for(path)
    for header, sequence in iter_entries(io.StringIO(text)):
        if not protein_name(header).startswith(ENTRAPMENT_PREFIX):
            continue
        for peptide in iter_peptides(sequence):
            if 7 <= len(peptide) <= 35:
                assert peptide not in originals, f"{peptide} collides with a real peptide"


def test_two_different_peptides_do_not_share_an_entrapment() -> None:
    """A collision between two distinct peptides must be retried, or one entrapment
    stands in for two originals and the pairing is ambiguous."""
    stats = EntrapmentStats()
    make_entrapment("ACDEFGHIKLMNPQRSTVWYK", stats=stats, originals=set())
    make_entrapment("STVWYACDEFGHIKLMNPQRK", stats=stats, originals=set())
    owners = list(stats._emitted.values())
    assert len(owners) == len(set(stats._emitted)), "an entrapment string served two peptides"


def test_a_repeated_peptide_is_not_treated_as_a_collision() -> None:
    """The subtle one. The rule is "no two DIFFERENT peptides share an entrapment", not
    "no entrapment string repeats". A naive already-emitted test would retry on the
    second occurrence of the same peptide and break the cross-database property."""
    stats = EntrapmentStats()
    peptide_protein = "ACDEFGHILMNPQTVWYK"
    first = make_entrapment(peptide_protein, stats=stats, originals=set())
    second = make_entrapment(peptide_protein, stats=stats, originals=set())
    assert first == second
    assert stats.distinct_peptides_retried == 0


def test_a_peptide_that_cannot_be_shuffled_is_counted_not_hidden() -> None:
    """A homopolymer has one arrangement. Keeping it silently would inflate N_E with
    peptides that are not verifiably false."""
    stats = EntrapmentStats()
    got = make_entrapment("AAAAAAAK", stats=stats, originals={"AAAAAAAK"})
    assert got == "AAAAAAAK"
    assert stats.peptides_unshuffled_collision == 1


def test_short_peptides_are_counted_separately_from_collisions() -> None:
    stats = EntrapmentStats()
    make_entrapment("AK", stats=stats, originals=set())
    assert stats.peptides_unshuffled_too_short == 1
    assert stats.peptides_unshuffled_collision == 0


def test_collect_peptides_keeps_only_the_verifiable_window(tmp_path: Path) -> None:
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": "AK" + "ACDEFGHIL" + "K"})
    peptides = collect_peptides(path)
    assert "AK" not in peptides
    assert all(7 <= len(p) <= 35 for p in peptides)


# --------------------------------------------------------------------------
# Ratio and pairing
# --------------------------------------------------------------------------
def test_ratio_one_writes_exactly_one_entrapment_per_protein(tmp_path: Path) -> None:
    entries = {f"sp|P{i:05d}|TEST": "ACDEFGHIKLMNPQRSTVWYK" for i in range(10)}
    path = write_fasta(tmp_path / "in.fasta", entries)
    names = [protein_name(h) for h, _ in iter_entries(io.StringIO(entrapments_for(path)))]
    assert sum(n.startswith(ENTRAPMENT_PREFIX) for n in names) == 10
    assert len(names) == 20


def test_a_fractional_ratio_selects_a_reproducible_subset(tmp_path: Path) -> None:
    entries = {f"sp|P{i:05d}|TEST": "ACDEFGHIKLMNPQRSTVWYK" for i in range(200)}
    path = write_fasta(tmp_path / "in.fasta", entries)
    first = entrapments_for(path, ratio=0.1)
    second = entrapments_for(path, ratio=0.1)
    assert first == second
    n = sum(
        protein_name(h).startswith(ENTRAPMENT_PREFIX) for h, _ in iter_entries(io.StringIO(first))
    )
    assert 5 <= n <= 35, f"0.1 of 200 proteins should be about 20, got {n}"


def test_ratio_zero_writes_the_database_unchanged(tmp_path: Path) -> None:
    entries = {"sp|P1|T": "ACDEFGHIKLMNPQRSTVWYK"}
    path = write_fasta(tmp_path / "in.fasta", entries)
    text = entrapments_for(path, ratio=0.0)
    assert ENTRAPMENT_PREFIX not in text


def test_the_entrapment_records_which_protein_it_came_from(tmp_path: Path) -> None:
    """ENTRAPMENT_DESIGN.md decision 7: nothing reads the pairing yet, but omitting it
    means regenerating the library, DIAMOND and all 12 runs to recover it."""
    path = write_fasta(tmp_path / "in.fasta", {"sp|P12345|FOO_HUMAN": "ACDEFGHIKLMNPQRK"})
    text = entrapments_for(path)
    assert f"{PAIR_TAG}sp|P12345|FOO_HUMAN" in text


def test_the_entrapment_accession_is_the_original_behind_the_prefix() -> None:
    """§2.3's prefix composition. GENERATE_LIBRARY_DECOYS applies the decoy prefix
    afterwards, giving LIBRARY_DECOY_ENTRAPMENT_sp|... -- so entrapment membership is
    tested by stripping the decoy prefix FIRST, never on a raw accession."""
    header = entrapment_header(">sp|P12345|FOO_HUMAN Some description", ENTRAPMENT_PREFIX, 0)
    assert protein_name(header) == "ENTRAPMENT_sp|P12345|FOO_HUMAN"
    decoyed = f"LIBRARY_DECOY_{protein_name(header)}"
    assert decoyed.removeprefix("LIBRARY_DECOY_").startswith(ENTRAPMENT_PREFIX)


def test_a_second_draw_gets_a_distinct_accession() -> None:
    first = entrapment_header(">sp|P1|T d", ENTRAPMENT_PREFIX, 0)
    second = entrapment_header(">sp|P1|T d", ENTRAPMENT_PREFIX, 1)
    assert protein_name(first) != protein_name(second)
    assert protein_name(second).startswith(ENTRAPMENT_PREFIX)


def test_the_original_description_is_kept(tmp_path: Path) -> None:
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T Cool protein OS=Mus": "ACDEFGHIKR"})
    assert "Cool protein OS=Mus" in entrapments_for(path)


# --------------------------------------------------------------------------
# Whole-file behaviour
# --------------------------------------------------------------------------
def test_every_original_entry_is_written_unchanged(tmp_path: Path) -> None:
    """Entrapments are added to the database, never substituted for it."""
    entries = {f"sp|P{i}|T": "ACDEFGHIKLMNPQRSTVWYK" for i in range(5)}
    path = write_fasta(tmp_path / "in.fasta", entries)
    out = {protein_name(h): s for h, s in iter_entries(io.StringIO(entrapments_for(path)))}
    for name, sequence in entries.items():
        assert out[name] == sequence


def test_a_trailing_stop_character_is_not_shuffled_into_the_sequence(tmp_path: Path) -> None:
    """Matches write_with_decoys. A '*' is not a residue and must not land mid-sequence."""
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": "ACDEFGHIKLMNPQRSTVWYK*"})
    for _header, sequence in iter_entries(io.StringIO(entrapments_for(path))):
        assert "*" not in sequence


def test_gzipped_input_round_trips(tmp_path: Path) -> None:
    """GENERATE_LIBRARY_DECOYS names its output .fasta.gz on this contract."""
    path = tmp_path / "in.fasta.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">sp|P1|T\nACDEFGHIKLMNPQRSTVWYK\n")
    with open_text(path) as handle:
        assert list(iter_entries(handle))
    text = entrapments_for(path)
    assert sum(1 for _ in iter_entries(io.StringIO(text))) == 2


def test_output_is_byte_identical_across_runs(tmp_path: Path) -> None:
    entries = {f"sp|P{i:05d}|TEST": "ACDEFGHIKLMNPQRSTVWYK" * 2 for i in range(30)}
    path = write_fasta(tmp_path / "in.fasta", entries)
    assert entrapments_for(path, seed=3) == entrapments_for(path, seed=3)


def test_a_different_seed_produces_a_different_database(tmp_path: Path) -> None:
    """Repeated draws are the escalation path for an inconclusive result, so two seeds
    must genuinely give two entrapment sets."""
    entries = {f"sp|P{i:05d}|TEST": "ACDEFGHIKLMNPQRSTVWYK" * 2 for i in range(30)}
    path = write_fasta(tmp_path / "in.fasta", entries)
    assert entrapments_for(path, seed=1) != entrapments_for(path, seed=2)


def test_stats_account_for_every_peptide(tmp_path: Path) -> None:
    """Counts must sum to the whole input, so 'not verifiably false' cannot hide."""
    entries = {f"sp|P{i}|T": "AK" + "ACDEFGHIKLMNPQRSTVWYK" * 2 for i in range(4)}
    path = write_fasta(tmp_path / "in.fasta", entries)
    buffer = io.StringIO()
    stats = write_with_entrapments(path, ENTRAPMENT_PREFIX, buffer, err=None)
    expected = sum(
        len(list(iter_peptides(s.rstrip("*"))))
        for _h, s in iter_entries(io.StringIO("".join(f">{k}\n{v}\n" for k, v in entries.items())))
    )
    counted = (
        stats.peptides_shuffled
        + stats.peptides_unshuffled_too_short
        + stats.peptides_unshuffled_collision
    )
    assert counted == expected
    assert stats.entrapments_written == 4


def test_the_stderr_report_names_the_verbatim_run_fraction(tmp_path: Path) -> None:
    """The headline diagnostic for the homology concern; it must not be silent."""
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": "ACDEFGHIKLMNPQRSTVWYK"})
    err = io.StringIO()
    write_with_entrapments(path, ENTRAPMENT_PREFIX, io.StringIO(), err=err)
    assert "verbatim run" in err.getvalue()


def test_a_collision_is_retried_with_a_different_shuffle() -> None:
    """The retry path. Seeding the shuffle on (peptide, seed, attempt) is what makes a
    retry explore instead of re-deriving the same colliding string forever."""
    peptide = "ACDEFGHILMNPQTVWYK"  # one tryptic peptide: no internal K or R
    first_try = shuffle_peptide(peptide, 0, 0)
    stats = EntrapmentStats()
    got = make_entrapment(peptide, stats=stats, originals={first_try})
    assert got != first_try
    assert got == shuffle_peptide(peptide, 0, 1)
    assert stats.distinct_peptides_retried == 1
    assert stats.peptides_shuffled == 1


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def test_cli_writes_originals_and_entrapments(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": "ACDEFGHIKLMNPQRSTVWYK"})
    assert generate_entrapments.main([str(path)]) == 0
    names = [protein_name(h) for h, _ in iter_entries(io.StringIO(capsys.readouterr().out))]
    assert names == ["sp|P1|T", "ENTRAPMENT_sp|P1|T"]


def test_cli_honours_the_entrapment_prefix(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": "ACDEFGHIKLMNPQRSTVWYK"})
    assert generate_entrapments.main([str(path), "--entrapment_prefix", "ENT2_"]) == 0
    assert "ENT2_sp|P1|T" in capsys.readouterr().out


def test_cli_gzipped_input_yields_gzipped_output(tmp_path: Path) -> None:
    """GENERATE_LIBRARY_DECOYS names its output .fasta.gz on this contract, so the
    entrapment step -- which runs immediately before it -- must honour it too."""
    path = tmp_path / "in.fasta.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">sp|P1|T\nACDEFGHIKLMNPQRSTVWYK\n")
    out = tmp_path / "out.fasta.gz"
    with out.open("wb") as sink:
        proc = subprocess.run(
            [sys.executable, "-m", "ms_denovo_db_utils.cli.generate_entrapments", str(path)],
            stdout=sink,
            stderr=subprocess.PIPE,
            check=True,
        )
    assert b"verbatim run" in proc.stderr
    with gzip.open(out, "rt") as handle:
        assert sum(1 for _ in iter_entries(handle)) == 2


def test_cli_can_reproduce_the_papers_literal_construction(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    long_peptide = "ACDEFGHILMNPQTVWY" * 3
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": long_peptide})
    assert generate_entrapments.main([str(path), "--outside_window", "keep"]) == 0
    sequences = [s for _h, s in iter_entries(io.StringIO(capsys.readouterr().out))]
    assert sequences == [long_peptide, long_peptide]


def test_cli_ratio_zero_adds_nothing(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": "ACDEFGHIKLMNPQRSTVWYK"})
    assert generate_entrapments.main([str(path), "--ratio", "0"]) == 0
    assert "ENTRAPMENT_" not in capsys.readouterr().out


# --------------------------------------------------------------------------
# Verbatim runs -- the homology diagnostic
# --------------------------------------------------------------------------
def test_a_verbatim_run_is_counted_only_when_it_is_long_enough() -> None:
    """Scattered fixed points are not homology. Only an uninterrupted stretch DIAMOND
    could seed on counts, which is why this is reported separately from residues_intact."""
    stats = EntrapmentStats()
    _record_verbatim_runs("ABCDEFGHIJ", "ABCDEFGHIJ", stats)
    assert stats.residues_in_verbatim_runs == 10
    assert stats.longest_verbatim_run == 10
    assert stats.proteins_with_verbatim_run == 1

    short = EntrapmentStats()
    _record_verbatim_runs("ABCDEFGHIJ", "ABCxxxxxxx", short)
    assert short.residues_in_verbatim_runs == 0
    assert short.longest_verbatim_run == 3
    assert short.proteins_with_verbatim_run == 0


def test_verbatim_run_totals_every_residue_of_a_long_run() -> None:
    stats = EntrapmentStats()
    _record_verbatim_runs("ABCDEFGHIJKL", "ABCDEFGHIJKL", stats)
    assert stats.residues_in_verbatim_runs == 12


def test_two_separate_runs_both_count() -> None:
    stats = EntrapmentStats()
    _record_verbatim_runs("ABCDEFGHxABCDEFGH", "ABCDEFGHyABCDEFGH", stats)
    assert stats.residues_in_verbatim_runs == 16
    assert stats.longest_verbatim_run == 8


def test_a_run_may_straddle_a_cleavage_boundary() -> None:
    """Measured over the whole protein, not per peptide. The residue before a K/R is the
    pinned C-terminus and the run continues past it, so per-peptide measurement would cut
    precisely the runs this is meant to find."""
    stats = EntrapmentStats()
    _record_verbatim_runs("ABCDKEFGHIJ", "ABCDKEFGHIJ", stats)
    assert stats.longest_verbatim_run == 11


def test_the_report_distinguishes_positional_matches_from_runs(tmp_path: Path) -> None:
    """The two numbers differ by an order of magnitude on real data and mean different
    things; the report must not let them be read as one."""
    path = write_fasta(tmp_path / "in.fasta", {"sp|P1|T": "ACDEFGHIKLMNPQRSTVWYK" * 3})
    err = io.StringIO()
    write_with_entrapments(path, ENTRAPMENT_PREFIX, io.StringIO(), err=err)
    text = err.getvalue()
    assert "hold their original position" in text
    assert "verbatim run" in text
