# ms-denovo-db-utils

Result-processing utilities for the [ms-denovo-db](https://github.com/mriffle/nf-ms-denovo-db)
Nextflow pipeline, and the container image that carries them.

The pipeline searches the same spectra with **Comet** (database search) and
**Casanovo** (de novo sequencing), homology-searches the union of the peptides
they report against a large annotated protein library with **DIAMOND**, and
estimates FDR over the resulting library identifications with
**percolator_RESET**. These tools are the glue between those steps.

## Tools

| Command | `/usr/local/bin` path | Role |
| --- | --- | --- |
| `process-comet-results` | `process_comet_results.py` | Collapse Comet PSMs to the best row per peptide, adding tryptic-terminus, ppm-error, decoy and rank features |
| `process-casanovo-results` | `process_casanovo_results.py` | The same for Casanovo mzTab, keyed on the highest score |
| `collate-into-fasta` | `collate_into_fasta.py` | Union both peptide lists into a DIAMOND query FASTA |
| `generate-reverse-decoys` | `generate_reverse_decoys.py` | Append a reversed decoy for every protein; gzip in, gzip out |
| `build-reset-input` | `build_reset_input.py` | Join Comet, Casanovo and DIAMOND results into the percolator_RESET feature table |

The pipeline invokes the absolute `/usr/local/bin/*.py` paths. Those are thin
shims in [`docker_bin/`](docker_bin/) that forward to `ms_denovo_db_utils.cli`;
the console-script names exist so the pipeline can migrate later.

### One row is one library peptide

`build-reset-input` inverts the usual unit of analysis. A row is not a PSM: it
is a *library peptide region*, the subject subsequence a DIAMOND alignment
landed on. Comet and Casanovo evidence for every query peptide aligning to that
region is aggregated onto it, and the target/decoy label comes from whether the
annotated-library protein carries the decoy prefix. FDR is therefore estimated
over annotated-database identifications, not spectrum identifications.

## Development

Requires Python 3.10+ and, for the container tests, Docker.

```sh
make install      # create .venv and install the package plus dev tools
make check        # everything CI runs
```

Individual targets: `lint`, `format`, `typecheck`, `test`, `image`,
`test-docker`. `make test-docker` skips with a message when no daemon is
available. Run `make help` for the full list.

CI runs these same targets, so a green `make check` locally means a green
build. Tests are held to `mypy --strict` alongside the package itself.

### Tests

- `tests/unit/` — behaviour of each module, plus end-to-end golden tests
- `tests/docker/` — the image's contract with Nextflow (script paths, entrypoint
  passthrough, exit-code propagation, `ps`, bash process substitution) and the
  same golden comparisons run through the built image
- `tests/fixtures/` — small hand-checkable inputs; see the README there

## Releasing

The version lives once, in `pyproject.toml`. Tagging `vX.Y.Z` runs a
build-and-verify workflow that refuses to proceed if the tag and the declared
version disagree. Publishing to quay.io is not wired up yet; see the note at
the bottom of `.github/workflows/release.yml`.
