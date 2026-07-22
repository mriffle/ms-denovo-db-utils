"""Guards on the packaging contract.

These are not smoke tests for their own sake. The release workflow refuses to
publish when the git tag disagrees with the declared version, and the container
image is built by ``pip install .`` -- so a version that cannot be read back
from installed metadata, or that drifts from pyproject.toml, breaks release
traceability (the problem OCI revision labels are meant to solve).
"""

from __future__ import annotations

import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - exercised only on the 3.10 CI leg
    import tomli as tomllib

import ms_denovo_db_utils

REPO_ROOT = Path(__file__).resolve().parents[2]


def _declared_version() -> str:
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        return str(tomllib.load(handle)["project"]["version"])


def test_installed_version_matches_pyproject() -> None:
    """The package must be installed, not merely importable off sys.path."""
    assert ms_denovo_db_utils.__version__ == _declared_version()


def test_version_is_not_the_uninstalled_fallback() -> None:
    assert ms_denovo_db_utils.__version__ != "0.0.0.dev0"


def test_package_supports_the_container_interpreter() -> None:
    """The image is python:3.10-slim; the package must not require anything newer."""
    with (REPO_ROOT / "pyproject.toml").open("rb") as handle:
        requires = str(tomllib.load(handle)["project"]["requires-python"])
    assert requires == ">=3.10"
    assert sys.version_info >= (3, 10)
