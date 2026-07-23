# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Five pure-standard-library Python tools that post-process search results for the
[nf-ms-denovo-db](https://github.com/mriffle/nf-ms-denovo-db) pipeline, plus the
container image (`quay.io/protio/ms-denovo-db-utils`) that ships them.

**Read [SPECIFICATION.md](SPECIFICATION.md) first.** It is the onboarding doc:
what each tool computes, column definitions, invariants, deliberate decisions
that look like bugs, and open TODOs. This file covers commands and the
essentials only.

## Commands

```sh
make install      # create .venv, install the package plus dev tools
make check        # everything CI runs: lint, typecheck, test, image, container tests
make format       # ruff check --fix + ruff format
```

Individual targets: `lint`, `typecheck`, `test`, `image`, `test-docker`, `clean`.
CI (`.github/workflows/ci.yml`) runs these same targets, so local and CI results
cannot drift.

Running specific tests — use the venv's pytest directly:

```sh
.venv/bin/pytest tests/unit/test_comet.py                       # one file
.venv/bin/pytest tests/unit/test_reset_input.py::test_name      # one test
.venv/bin/pytest -m "not docker"                                # skip container tests
.venv/bin/pytest -m docker                                      # container tests only
.venv/bin/pytest -m "not slow"                                  # skip the perf guard
```

`make test` adds `--cov` with a **95% floor**; a bare `pytest` invocation does not,
which is why running a single test does not fail on coverage.

## Essentials

- **No runtime dependencies, ever.** `dependencies = []` in `pyproject.toml` is
  deliberate — it is what keeps the image trivial to build and immune to
  dependency rot. Adding one needs a strong reason.
- **All logic lives in `src/ms_denovo_db_utils/`.** `docker_bin/` holds 3-line
  shims that exist solely to preserve the `/usr/local/bin/*.py` paths the
  pipeline invokes. Those paths and their **positional argument order are a
  deployed contract** — a pipeline revision may be pinned to an older image, so
  changes must be additive and optional, never reordered or removed.
- **Output must be deterministic.** Python randomises string hashing per process.
  Anything iterating a set (or a dict built from one) that reaches output must
  impose an order — **sort by default**. `tests/unit/test_determinism.py` runs the
  CLIs across eight `PYTHONHASHSEED` values. This has caused real defects twice;
  see SPECIFICATION.md §9.
- **A diff in `tests/fixtures/golden/` is a change in scientific output.** Those
  files were captured from the original flat scripts. Never regenerate one to make
  a test pass — work out why it changed, and if the change is correct update the
  golden in the same commit, with the reasoning in the commit message.
- **Malformed input must exit non-zero.** Nextflow reads exit codes; the original
  scripts printed a message and exited 0, producing an empty result file that read
  as success.
- **Python 3.10 is the contract** — it is what `python:3.10-slim` ships. The 3.11+
  CI legs are `continue-on-error` advisory warnings. `tomllib` is 3.11+, so code
  touching TOML needs the `tomli` fallback already in `tests/unit/test_packaging.py`.
- **`mypy --strict` covers `tests/` too** — new test functions need `-> None`.

## Footguns

- **A stale image gives a false green.** `pytest -m docker` reuses an existing
  `ms-denovo-db-utils:dev` and only builds when the tag is absent. After editing
  `src/`, run `make test-docker` (which rebuilds) or `make image` first.
- **`make check` without Docker is not full coverage.** `test-docker` prints `SKIP`
  and exits 0 by design, so a green run there says nothing about the image.
- **The image must keep `bash` and `procps`.** Every nf module writes through
  `> >(tee ...)` process substitution, which needs real bash; `ps` backs Nextflow's
  task metrics.
- **Gzip in means gzip out.** `generate_reverse_decoys` detects gzip by magic
  number, not extension, and `GENERATE_LIBRARY_DECOYS` in the pipeline names its
  output `.fasta.gz` or `.fasta` on exactly that behaviour.
- **The annotated library is streamed, never loaded whole.** It can be very large;
  `add_subject_sequences` reads one protein at a time on purpose.
