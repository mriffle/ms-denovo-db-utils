"""Generate entrapment proteins: shuffled sequences that cannot be in the sample.

Entrapment sequences are added to the searched database so the tool cannot tell them
from real ones. Any accepted hit on one is false by construction, and counting them
estimates the false discoveries hiding among the real targets. See
``ENTRAPMENT_DESIGN.md`` in the working directory, and Wen et al., *Nat Methods*
22(7):1454-1463, 2025.

Three things here are deliberate and easy to get wrong if edited:

**This must shuffle, never reverse.** :func:`~ms_denovo_db_utils.fasta.pseudo_reverse_sequence`
already implements the C-terminus-fixed, peptide-by-peptide idiom -- as *reversal*, to
build the library decoys. An entrapment generated the same way would be byte-identical to
the decoy of the same protein, and the experiment would be vacuous.

**Every peptide is shuffled, not just those in the length window.** This deviates from the
paper on purpose; see :func:`make_entrapment`.

**The same peptide must yield the same entrapment everywhere.** The shuffle is seeded from
the peptide sequence itself, so a peptide shared between the annotated library and the
Comet database gets the same entrapment in both -- the property
``pseudo_reverse_sequence``'s docstring explains for decoys. It uses ``hashlib`` rather
than ``hash()`` because ``hash()`` is salted per interpreter process and would silently
produce a different database on every run.
"""

from __future__ import annotations

import hashlib
import random
import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from .fasta import CLEAVAGE_RESIDUES, iter_entries, iter_peptides, open_text, protein_name

#: Peptide lengths the paper treats as identifiable. Outside this window a shuffled
#: peptide is not *verifiably* false, so it is excluded from the uniqueness check and
#: reported separately -- but it is still shuffled. See :func:`make_entrapment`.
MIN_VERIFIABLE_LENGTH = 7
MAX_VERIFIABLE_LENGTH = 35

#: How many times to re-shuffle a peptide that collided before giving up.
MAX_RETRIES = 20

#: Below this, a C-terminus-fixed shuffle has at most one movable position and cannot
#: produce anything new.
MIN_SHUFFLEABLE_LENGTH = 3

#: A contiguous verbatim stretch this long or longer is treated as an alignment anchor.
#: Set to the shortest peptide the pipeline considers identifiable: scattered single
#: residues that happen to stay put are not homology, an uninterrupted run is.
VERBATIM_RUN_THRESHOLD = 7

#: Marks the originating protein in an entrapment's description line, so the paired
#: estimator can be added later without regenerating the database (ENTRAPMENT_DESIGN.md
#: decision 7). Nothing reads it yet; that is the point.
PAIR_TAG = "entrapment_of="


#: Outcome categories for one peptide occurrence. Every occurrence lands in exactly one,
#: so the four counts sum to the number of peptides in the input.
SHUFFLED = "shuffled"
TOO_SHORT = "too_short"
COLLIDED = "collided"
KEPT_OUT_OF_WINDOW = "kept_out_of_window"


@dataclass
class EntrapmentStats:
    """What the generator could and could not guarantee, for the stderr report.

    The ``peptides_*`` counts are **per occurrence**, not per distinct sequence, and they
    partition the input: shuffled + too_short + collided + kept_out_of_window equals every
    tryptic peptide read. Occurrence is the right unit because these describe the database
    that was written -- a peptide left intact is intact everywhere it appears. The residue
    counts use the same unit, so the two can be read against each other.

    ``peptides_unshuffled_*`` are the numbers that matter scientifically: a peptide that
    could not be shuffled into something new is not verifiably false, and silently keeping
    it would inflate ``N_E``.

    **Two different "intact" numbers live here, and confusing them is easy.**
    ``residues_intact`` counts every position that happens to hold its original residue.
    It is large by construction and mostly uninformative: the C-terminus is pinned on
    purpose, so every shuffled peptide contributes at least one, and a random permutation
    contributes about one more. ``residues_in_verbatim_runs`` counts only positions inside
    an uninterrupted stretch of at least ``VERBATIM_RUN_THRESHOLD`` residues. **That is the
    one to read**, because a scattered fixed point is not something DIAMOND can align to
    and an uninterrupted run is.
    """

    proteins_read: int = 0
    entrapments_written: int = 0
    peptides_shuffled: int = 0
    peptides_unshuffled_too_short: int = 0
    peptides_unshuffled_collision: int = 0
    peptides_kept_out_of_window: int = 0
    peptides_outside_window: int = 0
    residues_total: int = 0
    residues_intact: int = 0
    residues_in_verbatim_runs: int = 0
    longest_verbatim_run: int = 0
    proteins_with_verbatim_run: int = 0
    #: Distinct peptides that needed at least one re-shuffle. A generation-time
    #: diagnostic, so distinct is the meaningful unit here -- unlike the counts above.
    distinct_peptides_retried: int = 0
    _cache: dict[str, tuple[str, str]] = field(default_factory=dict, repr=False)
    _emitted: dict[str, str] = field(default_factory=dict, repr=False)

    def count(self, category: str, peptide: str, entrapment: str) -> None:
        if category == SHUFFLED:
            self.peptides_shuffled += 1
        elif category == TOO_SHORT:
            self.peptides_unshuffled_too_short += 1
        elif category == COLLIDED:
            self.peptides_unshuffled_collision += 1
        else:
            self.peptides_kept_out_of_window += 1
        # strict=True: a shuffle is a permutation and the unshuffled paths return the
        # peptide itself, so a length change here would mean a real bug.
        self.residues_intact += sum(1 for a, b in zip(peptide, entrapment, strict=True) if a == b)

    def report(self, err: IO[str] = sys.stderr) -> None:
        total = self.residues_total
        intact = (100.0 * self.residues_intact / total) if total else 0.0
        print(
            f"entrapment: {self.entrapments_written} entrapment proteins from "
            f"{self.proteins_read} originals",
            file=err,
        )
        print(
            f"entrapment: {self.peptides_shuffled} peptide occurrences shuffled, "
            f"{self.distinct_peptides_retried} distinct peptides needed a retry",
            file=err,
        )
        print(
            f"entrapment: {self.peptides_unshuffled_too_short} peptides left intact "
            f"(<3 residues, nothing to permute), "
            f"{self.peptides_unshuffled_collision} left intact after failed retries, "
            f"{self.peptides_kept_out_of_window} kept verbatim as out-of-window "
            f"-- none of these are verifiably false",
            file=err,
        )
        print(
            f"entrapment: {self.peptides_outside_window} peptides outside the "
            f"{MIN_VERIFIABLE_LENGTH}-{MAX_VERIFIABLE_LENGTH} verifiable window "
            f"(excluded from the uniqueness check)",
            file=err,
        )
        anchored = (100.0 * self.residues_in_verbatim_runs / total) if total else 0.0
        print(
            f"entrapment: {self.residues_intact} / {self.residues_total} residues "
            f"({intact:.2f}%) hold their original position -- expected, the C-terminus "
            f"is pinned",
            file=err,
        )
        # The headline diagnostic. Entrapments are reached by homology in this pipeline,
        # so an uninterrupted verbatim stretch is what lets a real query peptide align to
        # an entrapment protein. Scattered fixed points, counted above, do not.
        print(
            f"entrapment: {self.residues_in_verbatim_runs} residues ({anchored:.2f}%) sit "
            f"in a verbatim run of >={VERBATIM_RUN_THRESHOLD} aa; "
            f"longest {self.longest_verbatim_run} aa; "
            f"{self.proteins_with_verbatim_run} / {self.entrapments_written} entrapments "
            f"carry one",
            file=err,
        )


def _peptide_rng(peptide: str, seed: int, attempt: int) -> random.Random:
    """A generator determined entirely by (peptide, seed, attempt).

    Not a stream RNG: the shuffle must not depend on which protein the peptide appeared
    in or where in the file it was, or the same peptide would get different entrapments
    in the annotated library and in the Comet database.
    """
    key = f"{seed}\x00{attempt}\x00{peptide}".encode()
    return random.Random(int.from_bytes(hashlib.sha256(key).digest()[:16], "big"))


def shuffle_peptide(peptide: str, seed: int = 0, attempt: int = 0) -> str:
    """Shuffle a peptide, holding its C-terminal residue fixed.

    Holding the C-terminus preserves the tryptic K/R, so the entrapment digests into
    peptides the same way the original does. A peptide of fewer than three residues has
    at most one movable position and is returned unchanged.
    """
    if len(peptide) < MIN_SHUFFLEABLE_LENGTH:
        return peptide
    body = list(peptide[:-1])
    _peptide_rng(peptide, seed, attempt).shuffle(body)
    return "".join(body) + peptide[-1]


def make_entrapment(  # noqa: PLR0913 - the paper's construction has this many knobs
    sequence: str,
    *,
    seed: int = 0,
    stats: EntrapmentStats | None = None,
    originals: frozenset[str] | set[str] | None = None,
    cleavage_residues: str = CLEAVAGE_RESIDUES,
    min_len: int = MIN_VERIFIABLE_LENGTH,
    max_len: int = MAX_VERIFIABLE_LENGTH,
    max_retries: int = MAX_RETRIES,
    shuffle_outside_window: bool = True,
) -> str:
    """Every tryptic peptide of ``sequence`` replaced by a shuffled counterpart.

    ``min_len``/``max_len`` bound the peptides the paper considers identifiable, and they
    govern the *uniqueness check*: only a peptide in the window is checked against the
    originals and against previously emitted entrapments, because only there does a
    collision matter for counting.

    **They do not govern whether a peptide is shuffled.** With
    ``shuffle_outside_window=True`` (the default) every peptide is shuffled, which is a
    deliberate deviation from the paper, for a reason specific to this pipeline:

    The paper's construction shuffles only in-window peptides, because in an ordinary
    database search an out-of-window peptide cannot be identified, so leaving it verbatim
    is harmless. **This pipeline does not match peptides exactly -- it reaches the library
    by homology (DIAMOND blastp).** Verbatim stretches are therefore not harmless: they
    are alignment anchors that let a *real* query peptide align to an entrapment protein,
    which either inflates ``N_E`` with true discoveries or dumps the region into the mixed
    bucket. Measured by this tool over Swiss-Prot minus mouse (556,421 proteins):

    =====================================  =========================  =====================
    ..                                     in-window only (`keep`)    every peptide
    =====================================  =========================  =====================
    residues in a verbatim run >= 7 aa     28.58%                     2.65%
    longest verbatim run                   5,048 aa                   168 aa
    entrapments carrying such a run        539,828 (97.0%)            313,324 (56.3%)
    =====================================  =========================  =====================

    Long tryptic peptides drive the difference: the window excludes peptides *above* 35
    residues as well as below 7, and those are exactly the long verbatim stretches DIAMOND
    aligns to. Pass ``shuffle_outside_window=False`` for the paper's literal construction.

    Note what the default does *not* fix. 56.3% of entrapments still carry some verbatim
    run of 7 aa or more, and the longest is 168 aa -- low-complexity regions, where every
    permutation is the original. Shuffling reduces the homology leak by an order of
    magnitude; it does not remove it, and §7.1's Stage B measurement is what bounds what
    is left.

    ``originals`` is the set of peptides a shuffle must avoid. Pass ``None`` to skip the
    check entirely (useful for tests and for single sequences).
    """
    stats = stats if stats is not None else EntrapmentStats()
    out: list[str] = []
    for peptide in iter_peptides(sequence, cleavage_residues):
        stats.residues_total += len(peptide)
        in_window = min_len <= len(peptide) <= max_len
        if not in_window:
            stats.peptides_outside_window += 1
            if not shuffle_outside_window:
                stats.count(KEPT_OUT_OF_WINDOW, peptide, peptide)
                out.append(peptide)
                continue
        out.append(
            _entrap_peptide(
                peptide,
                seed=seed,
                stats=stats,
                originals=originals,
                in_window=in_window,
                max_retries=max_retries,
            )
        )
    entrapment = "".join(out)
    # Measured over the whole protein, not per peptide: a verbatim run can straddle a
    # cleavage boundary, and the boundary residue is pinned on both sides of it, so
    # per-peptide measurement would cut exactly the runs it is meant to find.
    _record_verbatim_runs(sequence, entrapment, stats)
    return entrapment


def _record_verbatim_runs(original: str, entrapment: str, stats: EntrapmentStats) -> None:
    """Longest and total residues in uninterrupted verbatim stretches."""
    longest = run = total = 0
    for a, b in zip(original, entrapment, strict=True):
        run = run + 1 if a == b else 0
        longest = max(longest, run)
        if run >= VERBATIM_RUN_THRESHOLD:
            # The run only becomes reportable at the threshold; credit the residues that
            # got it there, then each further residue as it arrives.
            total += VERBATIM_RUN_THRESHOLD if run == VERBATIM_RUN_THRESHOLD else 1
    stats.residues_in_verbatim_runs += total
    stats.longest_verbatim_run = max(stats.longest_verbatim_run, longest)
    if longest >= VERBATIM_RUN_THRESHOLD:
        stats.proteins_with_verbatim_run += 1


def _entrap_peptide(  # noqa: PLR0913 - threads the caller's accumulators, not configuration
    peptide: str,
    *,
    seed: int,
    stats: EntrapmentStats,
    originals: frozenset[str] | set[str] | None,
    in_window: bool,
    max_retries: int,
) -> str:
    """One peptide's entrapment, with the uniqueness check and the repeat-peptide cache.

    The cache is not just an optimisation, it is what keeps the output well defined. The
    uniqueness rule is "no two *different* peptides may produce the same entrapment" --
    **not** "no entrapment string may repeat". The same peptide occurring again must
    produce the same entrapment, or the cross-database property this module exists to
    preserve is lost. A naive "have I emitted this string before?" test would retry on
    the second occurrence and break exactly that.
    """
    cached = stats._cache.get(peptide)
    if cached is not None:
        # Counted again: the counts are per occurrence, so that they partition the input
        # and can be read against the residue counts. Only the shuffling is memoised.
        entrapment, category = cached
        stats.count(category, peptide, entrapment)
        return entrapment

    def remember(entrapment: str, category: str) -> str:
        stats._cache[peptide] = (entrapment, category)
        stats.count(category, peptide, entrapment)
        return entrapment

    if len(peptide) < MIN_SHUFFLEABLE_LENGTH:
        return remember(peptide, TOO_SHORT)

    for attempt in range(max_retries):
        candidate = shuffle_peptide(peptide, seed, attempt)
        if not in_window or _acceptable(candidate, peptide, stats, originals):
            if attempt:
                stats.distinct_peptides_retried += 1
            if in_window:
                stats._emitted[candidate] = peptide
            return remember(candidate, SHUFFLED)

    # Could not find a new sequence -- short or low-complexity. Keep the original and
    # count it: it is not verifiably false, and hiding it would inflate N_E.
    return remember(peptide, COLLIDED)


def _acceptable(
    candidate: str,
    peptide: str,
    stats: EntrapmentStats,
    originals: frozenset[str] | set[str] | None,
) -> bool:
    """True if ``candidate`` collides with no real peptide and no other entrapment."""
    if originals is not None and candidate in originals:
        return False
    owner = stats._emitted.get(candidate)
    return owner is None or owner == peptide


def collect_peptides(
    input_file: str | Path,
    min_len: int = MIN_VERIFIABLE_LENGTH,
    max_len: int = MAX_VERIFIABLE_LENGTH,
    cleavage_residues: str = CLEAVAGE_RESIDUES,
) -> set[str]:
    """Every in-window tryptic peptide of the database, for the uniqueness check.

    A first pass over the input. Only in-window peptides are kept -- the check does not
    apply outside the window, and over Swiss-Prot this is the difference between 6.1M
    peptides and 22.9M.
    """
    peptides: set[str] = set()
    with open_text(input_file) as handle:
        for _header, sequence in iter_entries(handle):
            for peptide in iter_peptides(sequence.rstrip("*"), cleavage_residues):
                if min_len <= len(peptide) <= max_len:
                    peptides.add(peptide)
    return peptides


def _draws_for(name: str, ratio: float, seed: int) -> int:
    """How many entrapments this protein gets.

    The integer part applies to every protein; the fraction selects a reproducible random
    subset. Selection is keyed on the protein name rather than drawn from a stream, so it
    does not depend on file order and re-running yields the same subset.
    """
    whole = int(ratio)
    remainder = ratio - whole
    if remainder <= 0:
        return whole
    digest = hashlib.sha256(f"{seed}\x00select\x00{name}".encode()).digest()
    draw = int.from_bytes(digest[:8], "big") / float(1 << 64)
    return whole + (1 if draw < remainder else 0)


def entrapment_header(header: str, entrapment_prefix: str, index: int) -> str:
    """The entrapment's FASTA header, carrying its pairing back to the original.

    The accession is ``<prefix><original accession>``, so stripping the prefix recovers
    the original -- and, because the library decoy prefix is applied *afterwards* by
    ``GENERATE_LIBRARY_DECOYS``, a decoyed entrapment reads
    ``LIBRARY_DECOY_ENTRAPMENT_sp|...``. Test entrapment membership by stripping the decoy
    prefix first, never by testing ``ENTRAPMENT_`` on a raw accession.

    A second and subsequent draw appends ``_e<n>`` so the accessions stay unique without
    disturbing that prefix rule. The description repeats the pairing explicitly, because
    ``ratio > 1`` otherwise makes it recoverable only by convention.
    """
    original = protein_name(header)
    rest = header[1 + len(original) :].strip()
    suffix = "" if index == 0 else f"_e{index + 1}"
    description = f"{PAIR_TAG}{original}" + (f" {rest}" if rest else "")
    return f">{entrapment_prefix}{original}{suffix} {description}"


def iter_entrapments(  # noqa: PLR0913 - mirrors write_with_entrapments' signature
    input_file: str | Path,
    entrapment_prefix: str,
    originals: set[str] | None,
    *,
    ratio: float = 1.0,
    seed: int = 0,
    stats: EntrapmentStats | None = None,
    **kwargs: object,
) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` for every original entry, then its entrapments."""
    stats = stats if stats is not None else EntrapmentStats()
    with open_text(input_file) as handle:
        for header, sequence in iter_entries(handle):
            # A trailing stop character is not a residue and must not be shuffled into
            # the middle of the sequence. Matches write_with_decoys.
            trimmed = sequence.rstrip("*")
            stats.proteins_read += 1
            yield header, trimmed
            for index in range(_draws_for(protein_name(header), ratio, seed)):
                entrapment = make_entrapment(
                    trimmed,
                    seed=seed + index,
                    stats=stats,
                    originals=originals,
                    **kwargs,  # type: ignore[arg-type]
                )
                stats.entrapments_written += 1
                yield entrapment_header(header, entrapment_prefix, index), entrapment


def write_with_entrapments(  # noqa: PLR0913 - mirrors write_with_decoys plus the entrapment knobs
    input_file: str | Path,
    entrapment_prefix: str,
    output: IO[str],
    *,
    ratio: float = 1.0,
    seed: int = 0,
    err: IO[str] | None = sys.stderr,
    check_uniqueness: bool = True,
    **kwargs: object,
) -> EntrapmentStats:
    """Write every original entry, then ``ratio`` entrapment entries per original.

    Entrapments are injected **before** decoy generation, so that
    ``GENERATE_LIBRARY_DECOYS`` runs over the expanded database and every entrapment gets
    its own decoy. Injecting afterwards would leave targets outnumbering decoys and bias
    RESET's FDR low -- the one direction that makes the experiment useless.
    """
    originals = collect_peptides(input_file) if check_uniqueness else None
    stats = EntrapmentStats()
    for header, sequence in iter_entrapments(
        input_file, entrapment_prefix, originals, ratio=ratio, seed=seed, stats=stats, **kwargs
    ):
        output.write(f"{header}\n{sequence}\n")
    if err is not None:
        stats.report(err)
    return stats
