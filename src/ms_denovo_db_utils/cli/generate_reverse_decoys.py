"""Emit a FASTA containing every target protein followed by its reversed decoy."""

from __future__ import annotations

import argparse
import gzip
import sys
from collections.abc import Sequence

from ..fasta import is_gzipped, write_with_decoys


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Process FASTA file and generate decoy sequences.")
    parser.add_argument("input_file", help="Input FASTA file (can be gzipped)")
    parser.add_argument("--decoy_prefix", default="DECOY_", help="Decoy prefix (default: DECOY_)")
    args = parser.parse_args(argv)

    # A gzipped input yields a gzipped stdout; GENERATE_LIBRARY_DECOYS names its
    # output file .fasta.gz or .fasta on exactly this assumption.
    if is_gzipped(args.input_file):
        with gzip.open(sys.stdout.buffer, "wt") as output:
            write_with_decoys(args.input_file, args.decoy_prefix, output)
    else:
        write_with_decoys(args.input_file, args.decoy_prefix, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
