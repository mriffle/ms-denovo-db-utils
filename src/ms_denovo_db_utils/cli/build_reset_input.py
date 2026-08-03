"""Build the feature table consumed by percolator_RESET."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from ..diamond import add_subject_sequences, read_hits
from ..reset_input import read_casanovo_peptides, read_comet_peptides, write_reset_input


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build percolator_RESET input.")
    parser.add_argument("comet_results_file")
    parser.add_argument("casanovo_results_file")
    parser.add_argument("diamond_results_file")
    parser.add_argument("fasta_file", help="Annotated library FASTA, with decoys")
    parser.add_argument("library_decoy_prefix")
    # Unused: Comet decoy status is already resolved upstream into the is_decoy
    # column of comet_peptides.txt, so this prefix has nothing left to do here.
    #
    # Accepted and ignored permanently -- do NOT remove it. A pipeline revision
    # can be pinned to any published image, and `modules/build_reset_input.nf`
    # still passes this positionally, so dropping the parameter would break every
    # such run. Optional rather than required for the mirror case: a future
    # pipeline that stops passing it must still work against this image.
    #
    # (An earlier comment here promised removal "in the next commit". That was
    # never safe to do, for the reason above.)
    parser.add_argument("comet_decoy_prefix", nargs="?", default=None)
    # An OPTION, deliberately, not a seventh positional. The sixth positional above is
    # `nargs="?"`, so a seventh could never be distinguished from it being omitted --
    # `build_reset_input.py A B C D LIB ENT` would bind ENT to comet_decoy_prefix.
    # Empty means "no entrapment class", which is the pre-entrapment behaviour.
    parser.add_argument(
        "--entrapment_prefix",
        default=None,
        help="Accession prefix marking entrapment proteins in the annotated library "
        "(default: none, i.e. no entrapment accounting)",
    )
    args = parser.parse_args(argv)

    comet_map = read_comet_peptides(args.comet_results_file)
    casanovo_map = read_casanovo_peptides(args.casanovo_results_file)
    diamond_map = read_hits(
        args.diamond_results_file,
        args.library_decoy_prefix,
        entrapment_prefix=args.entrapment_prefix,
    )
    add_subject_sequences(diamond_map, args.fasta_file)

    write_reset_input(
        comet_map,
        casanovo_map,
        diamond_map,
        args.library_decoy_prefix,
        entrapment_prefix=args.entrapment_prefix,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
