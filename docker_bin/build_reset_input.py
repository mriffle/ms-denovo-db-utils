#!/usr/bin/env python3
"""Shim kept at /usr/local/bin/build_reset_input.py for the ms-denovo-db pipeline.

nf-ms-denovo-db invokes this path directly. The implementation lives in
ms_denovo_db_utils.cli.build_reset_input.
"""

import sys

from ms_denovo_db_utils.cli.build_reset_input import main

if __name__ == "__main__":
    sys.exit(main())
