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


def read_hits(file_path: str | Path, library_decoy_prefix: str) -> dict[str, DiamondHit]:
    """Return the best hit per query peptide, dropping ambiguous peptides.

    "Best" is the lowest e-value. A peptide that aligns to both target and
    decoy proteins in the annotated library carries no usable target/decoy
    signal, so it is discarded entirely -- regardless of which of its
    alignments happened to be the best one, or of the order they appear in.
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
        if _spans_target_and_decoy(proteins, library_decoy_prefix)
    }
    return {peptide: hit for peptide, hit in best.items() if peptide not in ambiguous}


def _spans_target_and_decoy(proteins: Iterable[str], library_decoy_prefix: str) -> bool:
    names = list(proteins)
    return any(n.startswith(library_decoy_prefix) for n in names) and any(
        not n.startswith(library_decoy_prefix) for n in names
    )


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
