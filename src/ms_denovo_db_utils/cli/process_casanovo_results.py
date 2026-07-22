"""Collapse Casanovo mzTab output to one best row per distinct peptide."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from ..casanovo import format_results, process_files


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Process Casanovo mzTab files.")
    parser.add_argument("files", nargs="+", help="Input Casanovo mzTab files")
    args = parser.parse_args(argv)

    peptides = process_files(args.files)
    for line in format_results(peptides):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
