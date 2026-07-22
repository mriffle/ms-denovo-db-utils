"""Collapse Comet PSM output to one best row per distinct peptide."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from ..comet import format_results, process_files


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Process Comet results files.")
    parser.add_argument(
        "--decoy_prefix",
        type=str,
        default="DECOY_",
        help="Decoy protein prefix (default: DECOY_)",
    )
    parser.add_argument("files", nargs="+", help="Input Comet result files")
    args = parser.parse_args(argv)

    peptides = process_files(args.files, args.decoy_prefix)
    for line in format_results(peptides):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
