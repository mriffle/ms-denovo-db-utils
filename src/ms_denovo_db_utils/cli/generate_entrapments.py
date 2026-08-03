"""Emit a FASTA containing every target protein followed by its entrapment(s).

Deliberately a separate command from ``generate_reverse_decoys``, rather than a flag on
it: the pipeline composes ``entrapments -> decoys`` as two steps, and keeping them
separate leaves that ordering visible in the DAG instead of buried in an option. The
order matters -- entrapments must exist before decoys are generated, so that every
entrapment gets its own decoy and target:decoy stays 1:1.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections.abc import Sequence

from ..entrapment import (
    MAX_RETRIES,
    MAX_VERIFIABLE_LENGTH,
    MIN_VERIFIABLE_LENGTH,
    write_with_entrapments,
)
from ..fasta import is_gzipped


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Add shuffled entrapment proteins to a FASTA, for FDR validation."
    )
    parser.add_argument("input_file", help="Input FASTA file (can be gzipped)")
    parser.add_argument(
        "--entrapment_prefix",
        default="ENTRAPMENT_",
        help="Prefix for entrapment accessions (default: ENTRAPMENT_)",
    )
    parser.add_argument(
        "--ratio",
        type=float,
        default=1.0,
        help=(
            "Entrapments per original protein (default: 1.0, i.e. paired). A fractional "
            "value adds an entrapment for that fraction of proteins, selected "
            "reproducibly from --seed. 1.0 is what the paired estimator needs."
        ),
    )
    parser.add_argument("--seed", type=int, default=0, help="Seed for the shuffles (default: 0)")
    parser.add_argument(
        "--min_len",
        type=int,
        default=MIN_VERIFIABLE_LENGTH,
        help=(
            f"Shortest peptide subject to the uniqueness check "
            f"(default: {MIN_VERIFIABLE_LENGTH}). Does NOT control what is shuffled."
        ),
    )
    parser.add_argument(
        "--max_len",
        type=int,
        default=MAX_VERIFIABLE_LENGTH,
        help=(
            f"Longest peptide subject to the uniqueness check "
            f"(default: {MAX_VERIFIABLE_LENGTH}). Does NOT control what is shuffled."
        ),
    )
    parser.add_argument(
        "--max_retries",
        type=int,
        default=MAX_RETRIES,
        help=f"Re-shuffles before giving up on a colliding peptide (default: {MAX_RETRIES})",
    )
    parser.add_argument(
        "--outside_window",
        choices=("shuffle", "keep"),
        default="shuffle",
        help=(
            "What to do with peptides outside --min_len..--max_len. 'shuffle' (default) "
            "shuffles them too; 'keep' is the paper's literal construction, which over "
            "Swiss-Prot leaves 28.6%% of residues in verbatim runs of >=7 aa (longest "
            "5,048 aa) against 2.65%% (longest 168 aa) when everything is shuffled. Those "
            "runs are harmless in an exact-match search and NOT harmless here, where the "
            "library is reached by homology."
        ),
    )
    parser.add_argument(
        "--no_uniqueness_check",
        action="store_true",
        help=(
            "Skip the first pass that collects every original peptide. Halves the run "
            "time and the memory, at the cost of allowing an entrapment peptide to "
            "coincide with a real one -- which would count a true identification as a "
            "false entrapment discovery."
        ),
    )
    args = parser.parse_args(argv)

    options = {
        "ratio": args.ratio,
        "seed": args.seed,
        "min_len": args.min_len,
        "max_len": args.max_len,
        "max_retries": args.max_retries,
        "shuffle_outside_window": args.outside_window == "shuffle",
        "check_uniqueness": not args.no_uniqueness_check,
    }

    # A gzipped input yields a gzipped stdout; GENERATE_LIBRARY_DECOYS names its output
    # file .fasta.gz or .fasta on exactly this assumption.
    if is_gzipped(args.input_file):
        with gzip.open(sys.stdout.buffer, "wt") as output:
            write_with_entrapments(args.input_file, args.entrapment_prefix, output, **options)
    else:
        write_with_entrapments(args.input_file, args.entrapment_prefix, sys.stdout, **options)
    return 0


if __name__ == "__main__":
    sys.exit(main())
