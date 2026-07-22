#!/usr/bin/env python3
"""Shim kept at /usr/local/bin/generate_reverse_decoys.py for the ms-denovo-db pipeline.

nf-ms-denovo-db invokes this path directly. The implementation lives in
ms_denovo_db_utils.cli.generate_reverse_decoys.
"""

import sys

from ms_denovo_db_utils.cli.generate_reverse_decoys import main

if __name__ == "__main__":
    sys.exit(main())
