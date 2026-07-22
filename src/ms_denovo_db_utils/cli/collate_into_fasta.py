"""Collect distinct peptides from result tables into a FASTA for homology search."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from ..collate import format_fasta, read_peptide_sequences


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Collate peptide sequences from result files into FASTA."
    )
    parser.add_argument("files", nargs="+", help="Tab-delimited files with peptides in column 1")
    args = parser.parse_args(argv)

    for line in format_fasta(read_peptide_sequences(args.files)):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
