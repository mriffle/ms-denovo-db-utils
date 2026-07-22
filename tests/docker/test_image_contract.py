"""The image's contract with the ms-denovo-db Nextflow pipeline.

Each assertion here corresponds to something nf-ms-denovo-db depends on and
would otherwise only discover during a pipeline run. These tests must keep
passing across the src/ restructuring -- that is what makes the move safe.
"""

from __future__ import annotations

import subprocess
from collections.abc import Callable

import pytest

pytestmark = pytest.mark.docker

Runner = Callable[..., "subprocess.CompletedProcess[str]"]

# nf-ms-denovo-db invokes these by absolute path:
#   modules/create_peptide_fasta.nf, modules/generate_decoys.nf,
#   modules/build_reset_input.nf
PIPELINE_SCRIPT_PATHS = [
    "/usr/local/bin/process_comet_results.py",
    "/usr/local/bin/process_casanovo_results.py",
    "/usr/local/bin/collate_into_fasta.py",
    "/usr/local/bin/build_reset_input.py",
    "/usr/local/bin/generate_reverse_decoys.py",
]


@pytest.mark.parametrize("script_path", PIPELINE_SCRIPT_PATHS)
def test_pipeline_script_path_exists(run_in_container: Runner, script_path: str) -> None:
    result = run_in_container(["test", "-f", script_path])
    assert result.returncode == 0, f"{script_path} is missing from the image"


@pytest.mark.parametrize("script_path", PIPELINE_SCRIPT_PATHS)
def test_pipeline_script_is_runnable(run_in_container: Runner, script_path: str) -> None:
    """Invoking with no arguments must fail cleanly, not crash on import.

    A missing module or syntax error surfaces here as a traceback rather than
    the scripts' own usage message.
    """
    result = run_in_container(["python3", script_path])
    combined = result.stdout + result.stderr
    assert "Traceback" not in combined, f"{script_path} raised on startup:\n{combined}"


def test_entrypoint_passes_the_command_through(run_in_container: Runner) -> None:
    """entrypoint.sh is `exec "$@"`; Nextflow relies on that passthrough."""
    result = run_in_container(["echo", "hello"])
    assert result.returncode == 0
    assert result.stdout.strip() == "hello"


def test_entrypoint_propagates_a_nonzero_exit_status(run_in_container: Runner) -> None:
    """Nextflow decides task success from the exit code, so it must survive."""
    result = run_in_container(["sh", "-c", "exit 42"])
    assert result.returncode == 42


def test_ps_is_available(run_in_container: Runner) -> None:
    """procps backs Nextflow's task resource metrics; losing it degrades reports."""
    result = run_in_container(["ps", "--version"])
    assert result.returncode == 0


def test_bash_process_substitution_works(run_in_container: Runner) -> None:
    """Every nf module writes output via `> >(tee ...)`, which needs real bash."""
    script = "echo payload > >(tee /tmp/captured.txt) ; sleep 0.2 ; cat /tmp/captured.txt"
    result = run_in_container(["bash", "-c", script])
    assert result.returncode == 0
    assert "payload" in result.stdout


def test_interpreter_is_python_310(run_in_container: Runner) -> None:
    """The base image pin; the package declares requires-python >=3.10."""
    result = run_in_container(["python3", "--version"])
    assert result.returncode == 0
    assert result.stdout.strip().startswith("Python 3.10.")
