#!/usr/bin/env python3
"""Shim kept at /usr/local/bin/process_casanovo_results.py for the ms-denovo-db pipeline.

nf-ms-denovo-db invokes this path directly. The implementation lives in
ms_denovo_db_utils.cli.process_casanovo_results.
"""

import sys

from ms_denovo_db_utils.cli.process_casanovo_results import main

if __name__ == "__main__":
    sys.exit(main())
