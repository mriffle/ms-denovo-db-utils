"""Result-processing utilities for the ms-denovo-db Nextflow pipeline."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("ms-denovo-db-utils")
except PackageNotFoundError:  # pragma: no cover - only when running from a raw checkout
    __version__ = "0.0.0.dev0"

__all__ = ["__version__"]
