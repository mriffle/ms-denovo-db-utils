# ms-denovo-db-utils

Result-processing utilities for the [ms-denovo-db](https://github.com/mriffle/nf-ms-denovo-db)
Nextflow pipeline, and the container image that carries them.

The pipeline searches the same mass spectra two ways — **Comet** (database
search) and **Casanovo** (de novo sequencing) — then homology-searches the
union of the peptides they report against a large annotated protein library
with **DIAMOND**, and estimates false discovery rate over the resulting library
identifications with **percolator_RESET**. The tools in this repository are the
glue between those steps: they collapse each search engine's output to
peptide-level evidence, build the homology-search query, generate decoys, and
assemble the feature table that drives the FDR estimate.

Everything here is pure Python standard library. There are no runtime
dependencies, which keeps the image quick to build and immune to dependency
rot.

## Purpose

Comet can only identify peptides present in the FASTA it was given, and
Casanovo can propose peptides that no database contains. Running both and
reconciling them against a large annotated library lets the pipeline ask a
different question than a normal search does: *which proteins in this large
annotated database are actually supported by the spectra?*

That reframing is why `build-reset-input` inverts the usual unit of analysis. A
row in its output is not a peptide-spectrum match — it is a **library peptide
region**, the subject subsequence a DIAMOND alignment landed on. Evidence from
every query peptide that aligned to that region is aggregated onto that one
row, and the target/decoy label comes from whether the annotated-library
protein carries the decoy prefix. FDR is therefore estimated over
annotated-database identifications rather than over spectrum identifications.

## What's in here

| Command | `/usr/local/bin` path | Role |
| --- | --- | --- |
| `process-comet-results` | `process_comet_results.py` | Collapse Comet PSMs to the best row per peptide, adding tryptic-terminus, ppm-error, decoy and rank features |
| `process-casanovo-results` | `process_casanovo_results.py` | The same for Casanovo mzTab, keyed on the highest score |
| `collate-into-fasta` | `collate_into_fasta.py` | Union both peptide lists into a DIAMOND query FASTA |
| `generate-reverse-decoys` | `generate_reverse_decoys.py` | Append a reversed decoy for every protein; gzip in, gzip out |
| `build-reset-input` | `build_reset_input.py` | Join Comet, Casanovo and DIAMOND results into the percolator_RESET feature table |

```
src/ms_denovo_db_utils/   the package — all logic lives here
docker_bin/               thin shims installed to /usr/local/bin/*.py
tests/unit/               module behaviour, golden, determinism, performance
tests/docker/             image contract, and golden output through the image
tests/fixtures/           small hand-checkable inputs and expected outputs
Dockerfile, entrypoint.sh the image
Makefile                  the development harness; CI runs the same targets
```

The Nextflow pipeline invokes the absolute `/usr/local/bin/*.py` paths. Those
are shims in [`docker_bin/`](docker_bin/) that forward into
`ms_denovo_db_utils.cli`; the console-script names in the table above exist so
the pipeline can migrate to plain commands whenever convenient.

For the full picture — data formats, column definitions, invariants, known
issues and outstanding work — see [`SPECIFICATION.md`](SPECIFICATION.md).

## Building the Docker image

The image is `python:3.10-slim` plus `procps` (Nextflow shells out to `ps` for
task resource metrics), with the package `pip install`ed and the shims placed
on `/usr/local/bin`.

```sh
docker build -t ms-denovo-db-utils:dev .
```

or equivalently `make image`.

### Build arguments

Both are optional and default to placeholders. They populate OCI labels, which
is what lets a running container be traced back to a commit.

| Argument | Meaning |
| --- | --- |
| `VERSION` | Image version, normally the git tag |
| `VCS_REF` | Commit SHA the image was built from |

A release-style build:

```sh
docker build \
  --build-arg VERSION=1.2.0 \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  -t quay.io/protio/ms-denovo-db-utils:1.2.0 .
```

Inspect the result:

```sh
docker inspect --format '{{ json .Config.Labels }}' ms-denovo-db-utils:dev
```

### Running a tool from the image

The entrypoint is a passthrough (`exec "$@"`), so any command works and exit
codes propagate — which is what Nextflow reads to decide task success.

```sh
docker run --rm -v "$PWD":/data -w /data ms-denovo-db-utils:dev \
  python3 /usr/local/bin/process_comet_results.py \
    --decoy_prefix COMET_DECOY_ run1.txt run2.txt > comet_peptides.txt
```

Note that the pipeline writes output through bash process substitution
(`> >(tee ...)`), so the image must keep a real `bash`. There is a test for it.

### Testing the image

```sh
make test-docker
```

This rebuilds the image and then runs the container suite, which checks the
contract the pipeline depends on — that all five script paths exist and start
cleanly, that the entrypoint passes commands through and propagates non-zero
exit codes, that `ps` and bash process substitution work — and re-runs the full
tool chain inside the container, comparing against the same expected outputs
the unit tests assert on. That last part is what catches a broken install layer
or a shim that silently fails.

> After changing anything under `src/`, rebuild before running container tests.
> `pytest -m docker` on its own reuses an existing `ms-denovo-db-utils:dev`
> image and will happily give you a green run against stale code. `make
> test-docker` rebuilds for you. Point the suite at a different tag with the
> `MS_DENOVO_DB_UTILS_IMAGE` environment variable.

### Publishing

Not wired up yet: this repository is not linked to quay.io. Tagging `vX.Y.Z`
runs a build-and-verify workflow that refuses to proceed when the git tag and
the version declared in `pyproject.toml` disagree. The login and push steps to
add once a robot account exists are written out at the bottom of
[`.github/workflows/release.yml`](.github/workflows/release.yml).

The pipeline consumes whatever tag is named by
`params.images.ms_denovo_db_utils` in `nf-ms-denovo-db/container_images.config`,
so publishing a new image is only half of a release — the pipeline has to be
pointed at it.

## Development

Requires Python 3.10+ and, for the container tests, Docker.

```sh
make install      # create .venv and install the package plus dev tools
make check        # everything CI runs
```

Individual targets: `lint`, `format`, `typecheck`, `test`, `image`,
`test-docker`, `clean`. Run `make help` for the list. `make test-docker` skips
with a message when no daemon is available, so `make check` stays usable on a
machine without Docker — which also means a green run there says nothing about
the image.

CI runs these same targets, so a green `make check` locally means a green
build. The package and the tests are both held to `mypy --strict`.

Tests are described in [`SPECIFICATION.md`](SPECIFICATION.md#8-development).
One rule worth stating here: the files in `tests/fixtures/golden/` are expected
scientific output, captured before this repository was restructured. A diff in
one of them means the pipeline now produces something different. Never
regenerate a golden file to make a test pass.

## License

Apache License 2.0 — see [LICENSE](LICENSE).
