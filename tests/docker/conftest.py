"""Fixtures for tests that exercise the built container image.

Every test in this package is marked ``docker`` and is skipped -- loudly, not
silently -- when no usable daemon is present, so ``make check`` still works on
a machine without Docker.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
IMAGE_ENV_VAR = "MS_DENOVO_DB_UTILS_IMAGE"
DEFAULT_IMAGE = "ms-denovo-db-utils:dev"


def _daemon_is_usable() -> bool:
    if shutil.which("docker") is None:
        return False
    return subprocess.run(["docker", "info"], capture_output=True, check=False).returncode == 0


def _image_exists(tag: str) -> bool:
    result = subprocess.run(["docker", "image", "inspect", tag], capture_output=True, check=False)
    return result.returncode == 0


@pytest.fixture(scope="session")
def container_image() -> str:
    """Tag of the image under test, building it if it is not already present.

    CI builds the image with a cached buildx step before invoking pytest; this
    fixture only pays the build cost on a developer machine.
    """
    if not _daemon_is_usable():
        pytest.skip("no usable Docker daemon")

    tag = os.environ.get(IMAGE_ENV_VAR, DEFAULT_IMAGE)
    if not _image_exists(tag):
        build = subprocess.run(
            ["docker", "build", "-t", tag, str(REPO_ROOT)],
            capture_output=True,
            text=True,
            check=False,
        )
        if build.returncode != 0:
            pytest.fail(f"docker build failed:\n{build.stdout}\n{build.stderr}")
    return tag


@pytest.fixture(scope="session")
def run_in_container(
    container_image: str,
) -> Callable[..., subprocess.CompletedProcess[str]]:
    """Return a helper that runs a command inside the image under test."""

    def _run(
        argv: Sequence[str],
        *,
        mounts: dict[Path, str] | None = None,
        workdir: str | None = None,
        as_current_user: bool = False,
        env: dict[str, str] | None = None,
        check: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        command = ["docker", "run", "--rm"]
        for host_path, container_path in (mounts or {}).items():
            command += ["-v", f"{host_path.resolve()}:{container_path}"]
        if as_current_user:
            # Anything written into a mounted directory would otherwise be
            # owned by root and defeat pytest's tmp_path cleanup.
            command += ["-u", f"{os.getuid()}:{os.getgid()}"]
        for name, value in (env or {}).items():
            command += ["-e", f"{name}={value}"]
        if workdir is not None:
            command += ["-w", workdir]
        command.append(container_image)
        command += list(argv)
        return subprocess.run(command, capture_output=True, text=True, check=check)

    return _run
