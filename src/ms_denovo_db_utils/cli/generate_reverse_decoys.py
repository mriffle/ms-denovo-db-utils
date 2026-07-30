"""Emit a FASTA containing every target protein followed by its reversed decoy."""

from __future__ import annotations

import argparse
import gzip
import sys
from collections.abc import Sequence

from ..fasta import PEPTIDE_REVERSAL, REVERSAL_METHODS, is_gzipped, write_with_decoys


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Process FASTA file and generate decoy sequences.")
    parser.add_argument("input_file", help="Input FASTA file (can be gzipped)")
    parser.add_argument("--decoy_prefix", default="DECOY_", help="Decoy prefix (default: DECOY_)")
    parser.add_argument(
        "--reversal_method",
        choices=REVERSAL_METHODS,
        default=PEPTIDE_REVERSAL,
        help=(
            "How to derive a decoy. 'peptide' (default) reverses each tryptic "
            "peptide with its C-terminal residue held in place, matching how the "
            "custom Casanovo model's training decoys were built. 'protein' "
            "reverses the whole sequence, which is what this tool did before."
        ),
    )
    args = parser.parse_args(argv)

    # A gzipped input yields a gzipped stdout; GENERATE_LIBRARY_DECOYS names its
    # output file .fasta.gz or .fasta on exactly this assumption.
    if is_gzipped(args.input_file):
        with gzip.open(sys.stdout.buffer, "wt") as output:
            write_with_decoys(args.input_file, args.decoy_prefix, output, args.reversal_method)
    else:
        write_with_decoys(args.input_file, args.decoy_prefix, sys.stdout, args.reversal_method)
    return 0


if __name__ == "__main__":
    sys.exit(main())
