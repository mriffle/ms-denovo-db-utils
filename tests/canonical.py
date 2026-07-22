"""Normalisation of environment-dependent noise in tool output.

Result tables record the input path they were built from, which differs between
a local run and a run inside the container. Only the basename is meaningful, so
that column is normalised before comparing against a golden file.

Row *ordering* used to need normalising too, because output came out in Python
set-iteration order. It no longer does: the tools emit a defined order, and the
tests assert raw equality.
"""

from __future__ import annotations

from pathlib import PurePosixPath


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
