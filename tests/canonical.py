"""Normalisers that make tool output comparable to a stored golden file.

Two kinds of noise are removed:

*Environment noise* -- result tables record the input path they came from, which
differs between a local run and a run inside the container. Only the basename
is meaningful.

*Ordering noise* -- the tools currently emit rows in Python set-iteration order,
so row order and ScanNr vary between runs. The canonical order defined here
(rows sorted by SpecId, ScanNr renumbered to match) is the order the tools are
made to emit natively once ordering is fixed, at which point these helpers
become identity transforms on well-formed output rather than repairs.
"""

from __future__ import annotations

from pathlib import PurePosixPath

SPEC_ID_COLUMN = 0
SCAN_NR_COLUMN = 2


def basename_column(text: str, column: int) -> str:
    """Replace a path-valued column with its basename in every data row."""
    lines = text.splitlines()
    if not lines:
        return text

    out = [lines[0]]
    for line in lines[1:]:
        fields = line.split("\t")
        if column < len(fields):
            fields[column] = PurePosixPath(fields[column]).name
        out.append("\t".join(fields))
    return "\n".join(out) + "\n"


def canonical_reset_input(text: str) -> str:
    """Sort reset_input rows by SpecId and renumber ScanNr to match."""
    lines = [line for line in text.splitlines() if line]
    if not lines:
        return text

    header, rows = lines[0], lines[1:]
    parsed = sorted((row.split("\t") for row in rows), key=lambda f: f[SPEC_ID_COLUMN])
    for scan_nr, fields in enumerate(parsed, start=1):
        fields[SCAN_NR_COLUMN] = str(scan_nr)
    return "\n".join([header, *("\t".join(f) for f in parsed)]) + "\n"


def canonical_fasta(text: str) -> str:
    """Sort FASTA entries by header so set-ordered output compares equal."""
    lines = [line for line in text.splitlines() if line]
    entries = sorted((lines[i], lines[i + 1]) for i in range(0, len(lines) - 1, 2))
    return "\n".join(f"{header}\n{sequence}" for header, sequence in entries) + "\n"


def canonical_missing_warning(stderr: str) -> str:
    """Sort the peptide list in the 'missing from diamond_map' warning."""
    marker = "will be skipped: "
    out = []
    for line in stderr.splitlines():
        if marker in line:
            prefix, _, peptides = line.partition(marker)
            listed = sorted(p.strip() for p in peptides.split(","))
            out.append(prefix + marker + ", ".join(listed))
        else:
            out.append(line)
    return "\n".join(out) + "\n" if out else ""
