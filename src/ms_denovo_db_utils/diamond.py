"""Read DIAMOND ``--outfmt 6`` results and resolve the subject sequences hit."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from .fasta import iter_entries, open_text, protein_name

#: DIAMOND's tabular default, in order.
COLUMN_HEADERS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)


@dataclass
class DiamondHit:
    """The best alignment found for one query peptide."""

    qseqid: str
    sseqid: str
    pident: float
    sstart: int
    send: int
    evalue: float
    bitscore: float
    #: Residues of the subject protein covered by the alignment. Filled in by
    #: :func:`add_subject_sequences`; this is the "library peptide" that
    #: reset_input groups on.
    ssequence: str | None = None

    @property
    def subject_length(self) -> int:
        return self.send - self.sstart + 1

    def is_decoy(self, library_decoy_prefix: str) -> bool:
        return self.sseqid.startswith(library_decoy_prefix)


def parse_line(line: str, line_num: int) -> DiamondHit:
    columns = line.split("\t")
    if len(columns) != len(COLUMN_HEADERS):
        raise ValueError(
            f"Line {line_num}: Expected {len(COLUMN_HEADERS)} columns, but found {len(columns)}"
        )
    return DiamondHit(
        qseqid=columns[0],
        sseqid=columns[1],
        pident=float(columns[2]),
        sstart=int(columns[8]),
        send=int(columns[9]),
        evalue=float(columns[10]),
        bitscore=float(columns[11]),
    )


def _preference(hit: DiamondHit) -> tuple[float, float, str]:
    """Sort key selecting the best hit; smaller is better.

    Ranked by e-value, then by descending bit score, then by subject name so
    that fully tied alignments resolve the same way no matter what order
    DIAMOND emitted its rows in.
    """
    return (hit.evalue, -hit.bitscore, hit.sseqid)


#: The four classes a library subject can belong to. `Label` is decided by the
#: decoy axis alone; entrapment membership is orthogonal to it.
TARGET = "target"
ENTRAPMENT_TARGET = "entrapment_target"
DECOY = "decoy"
ENTRAPMENT_DECOY = "entrapment_decoy"


def classify_subject(
    sseqid: str, library_decoy_prefix: str, entrapment_prefix: str | None = None
) -> str:
    """Which of the four classes a library subject belongs to.

    **Strip the decoy prefix first, then test for entrapment.** The prefixes compose with
    the decoy prefix outermost, because entrapments are injected *before* decoy generation
    and so the decoy step wraps them::

        sp|P12345|FOO                           -> target
        ENTRAPMENT_sp|P12345|FOO                -> entrapment_target
        LIBRARY_DECOY_sp|P12345|FOO             -> decoy
        LIBRARY_DECOY_ENTRAPMENT_sp|P12345|FOO  -> entrapment_decoy

    Testing ``ENTRAPMENT_`` against a raw accession would classify every decoyed
    entrapment as a plain decoy.
    """
    is_decoy = sseqid.startswith(library_decoy_prefix)
    remainder = sseqid[len(library_decoy_prefix) :] if is_decoy else sseqid
    is_entrapment = bool(entrapment_prefix) and remainder.startswith(str(entrapment_prefix))
    if is_decoy:
        return ENTRAPMENT_DECOY if is_entrapment else DECOY
    return ENTRAPMENT_TARGET if is_entrapment else TARGET


def read_hits(
    file_path: str | Path,
    library_decoy_prefix: str,
    *,
    entrapment_prefix: str | None = None,
) -> dict[str, DiamondHit]:
    """Return the best hit per query peptide, dropping ambiguous peptides.

    "Best" is the lowest e-value. A peptide that aligns to both target and
    decoy proteins in the annotated library carries no usable target/decoy
    signal, so it is discarded entirely -- regardless of which of its
    alignments happened to be the best one, or of the order they appear in.

    ``entrapment_prefix`` does **not** change which peptides are dropped: the drop rule is
    the target/decoy axis, and entrapments are targets. It exists so that a caller can ask
    for the entrapment classes; with it ``None`` this function is bit-identical to what it
    was before entrapment existed.
    """
    best: dict[str, DiamondHit] = {}
    proteins_by_peptide: dict[str, set[str]] = {}

    with Path(file_path).open() as handle:
        for line_num, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            hit = parse_line(line, line_num)

            incumbent = best.get(hit.qseqid)
            if incumbent is None or _preference(hit) < _preference(incumbent):
                best[hit.qseqid] = hit

            # Accumulated over *every* alignment, so the ambiguity check below
            # sees the peptide's full protein set.
            proteins_by_peptide.setdefault(hit.qseqid, set()).add(hit.sseqid)

    ambiguous = {
        peptide
        for peptide, proteins in proteins_by_peptide.items()
        if _spans_target_and_decoy(proteins, library_decoy_prefix, entrapment_prefix)
    }
    return {peptide: hit for peptide, hit in best.items() if peptide not in ambiguous}


def _spans_target_and_decoy(
    proteins: Iterable[str], library_decoy_prefix: str, entrapment_prefix: str | None = None
) -> bool:
    """True iff a peptide's alignments span the target/decoy axis.

    Written through :func:`classify_subject` so the drop rule is explicit in terms of the
    four classes rather than a bare prefix test that silently collapses them:
    ``entrapment_target`` counts as a target and ``entrapment_decoy`` as a decoy.

    **Deliberately blind to the entrapment axis.** A query reaching both an original and
    an entrapment protein of the *same* label is not dropped -- both are targets, and
    dropping them would remove exactly the evidence the experiment exists to count. Such a
    region is marked *mixed* by :func:`region_entrapment_class` instead, and excluded from
    N_T and N_E there, because it is not cleanly attributable to either population.

    Passing ``entrapment_prefix`` therefore cannot change the result: it only ever splits a
    class into two classes that map to the same side of this test.
    """
    decoy_flags = {
        classify_subject(name, library_decoy_prefix, entrapment_prefix) in (DECOY, ENTRAPMENT_DECOY)
        for name in proteins
    }
    return len(decoy_flags) > 1


#: Positions on the entrapment axis. `MIXED` is a region reachable from both an original
#: and an entrapment protein: not attributable to either population, so it is counted
#: separately and excluded from N_T and N_E.
ORIGINAL = "original"
ENTRAPMENT = "entrapment"
MIXED = "mixed"


def region_entrapment_axis(
    sseqids: Iterable[str], library_decoy_prefix: str, entrapment_prefix: str | None = None
) -> str:
    """Where a library *region* sits on the entrapment axis: original, entrapment or mixed.

    A row of reset_input.txt is a region, and its `Proteins` column names the whole
    retained group -- so a region can be reachable from an original protein and from an
    entrapment at once. That region cannot be counted as an entrapment discovery without
    also counting it as a target one, so it gets its own class.

    **This reports the entrapment axis only, never the target/decoy one.** A region can
    name both a target and a decoy protein -- two different query peptides, each
    individually unambiguous, landing on the same subsequence in a target and in a decoy
    protein. That is pre-existing, documented (SPECIFICATION.md §5), and has nothing to do
    with entrapment: measured over the mouse benchmark it is 2 regions of 196,085, on a
    library containing no entrapments at all. Folding it in here would pollute the one
    count that decides whether the entrapment experiment is compromised. The row's
    target/decoy side comes from its best hit, which is what `Label` already means.
    """
    flags = {
        classify_subject(name, library_decoy_prefix, entrapment_prefix)
        in (ENTRAPMENT_TARGET, ENTRAPMENT_DECOY)
        for name in sseqids
    }
    if not flags:
        raise ValueError("a region must name at least one protein")
    if len(flags) > 1:
        return MIXED
    return ENTRAPMENT if next(iter(flags)) else ORIGINAL


def region_entrapment_class(
    sseqids: Iterable[str],
    library_decoy_prefix: str,
    entrapment_prefix: str | None = None,
    *,
    is_decoy_row: bool = False,
) -> str:
    """The reporting class of a region: one of the four, or ``mixed``.

    ``is_decoy_row`` is the region's `Label`, which is decided by its best hit alone --
    not re-derived from the protein list, so this cannot disagree with the row.
    """
    axis = region_entrapment_axis(sseqids, library_decoy_prefix, entrapment_prefix)
    if axis == MIXED:
        return MIXED
    if is_decoy_row:
        return ENTRAPMENT_DECOY if axis == ENTRAPMENT else DECOY
    return ENTRAPMENT_TARGET if axis == ENTRAPMENT else TARGET


def add_subject_sequences(hits: dict[str, DiamondHit], fasta_file_path: str | Path) -> None:
    """Fill in :attr:`DiamondHit.ssequence` for every hit, in place.

    The annotated library can be very large, so it is streamed one protein at a
    time rather than loaded into memory.
    """
    hits_by_protein: dict[str, list[DiamondHit]] = {}
    for hit in hits.values():
        hits_by_protein.setdefault(hit.sseqid, []).append(hit)

    with open_text(fasta_file_path) as handle:
        for header, sequence in iter_entries(handle):
            for hit in hits_by_protein.get(protein_name(header), ()):
                hit.ssequence = sequence[hit.sstart - 1 : hit.send]

    missing = [peptide for peptide, hit in hits.items() if hit.ssequence is None]
    if missing:
        raise ValueError(f"Peptide {missing[0]} does not have an 'ssequence' property")
