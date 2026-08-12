"""Peptide identity, where two spellings mean one molecule."""

from __future__ import annotations


def normalise_il(sequence: str) -> str:
    """Canonical spelling of a peptide under I == L.

    Leucine and isoleucine are structural isomers of identical residue mass, so no
    mass spectrometer distinguishes them and neither engine can. Casanovo is run with
    ``replace_isoleucine_with_leucine: true`` and emits **zero** I; Comet reports
    whichever letter the database spells, and 49% of its peptides contain an I. The
    two engines therefore cannot spell the same observed peptide the same way whenever
    an I is involved, and without this the I-form and the L-form are two separate
    hypotheses in the FDR pool.

    **L is the canonical form because it is the one Casanovo already emits.** Choosing
    I would mean rewriting every de novo peptide to agree with the database ones rather
    than the other way round, and the de novo side is the one with no information to
    lose -- its L is already a convention, not a measurement.

    This is applied when the two engines' results are **merged**, not when the query
    FASTA is written. Both spellings stay DIAMOND queries so that each aligns against
    the un-normalised library at its own best spelling: normalising the query alone
    would turn 41,342 exact I<->I matches into I<->L substitutions across 24,691 Comet
    peptides that currently align at 100% identity, degrading two RESET features for
    roughly half the database peptides in order to merge 1,606.
    """
    return sequence.replace("I", "L")
