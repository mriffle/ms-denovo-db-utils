# SPECIFICATION

Onboarding document for anyone — human or agent — picking up work on this
repository. It covers what the project is, how it works, how it is developed,
what is deliberately the way it is, and what is left to do.

Read [`README.md`](README.md) first for the user-facing overview. This document
is the operational detail behind it.

---

## 1. What this is

`ms-denovo-db-utils` is the source repository for the container image
`quay.io/protio/ms-denovo-db-utils`, which supplies five result-processing
tools to the [ms-denovo-db](https://github.com/mriffle/nf-ms-denovo-db)
Nextflow pipeline.

The tools are pure Python standard library — **there are no runtime
dependencies and none should be added without a strong reason.** That is what
keeps the image trivial to build and immune to dependency rot.

### Relationship to the pipeline repository

Two separate repositories, coupled by an image tag:

| | |
| --- | --- |
| This repo | `github.com/mriffle/ms-denovo-db-utils` — the tools and the image |
| Pipeline | `github.com/mriffle/nf-ms-denovo-db` — the Nextflow workflow |
| Coupling | `params.images.ms_denovo_db_utils` in `nf-ms-denovo-db/container_images.config` |

On a development machine the pipeline is typically checked out alongside this
repo at `../nf-ms-denovo-db`.

**The pipeline calls these tools by absolute path**, from three modules:

| nf module | Calls |
| --- | --- |
| `modules/create_peptide_fasta.nf` | `process_comet_results.py`, `process_casanovo_results.py`, `collate_into_fasta.py` |
| `modules/generate_decoys.nf` | `generate_reverse_decoys.py` (twice: Comet DB and annotated library) |
| `modules/build_reset_input.nf` | `build_reset_input.py` |

Those `/usr/local/bin/*.py` paths are a hard interface. Breaking them breaks a
pipeline that may be pinned to an older image tag and cannot be fixed in
lockstep. See §7.

---

## 2. Current status

| | |
| --- | --- |
| Branch | `main` (pushed; `test-harness-and-fixes` is a leftover local branch, identical) |
| Version | `1.2.0.dev0` in `pyproject.toml` — **not yet released** |
| Image the pipeline currently pins | `quay.io/protio/ms-denovo-db-utils:1.1.10` (pre-dates all of this work) |
| Tests | 190 unit + 22 container |
| Coverage | ~99%, floor enforced at 95% |
| Type checking | `mypy --strict` over `src/` **and** `tests/`, clean |
| Lint | `ruff` on a broad rule set, clean; `hadolint` on the Dockerfile, clean |

### What has been done

1. Development harness: `pyproject.toml`, `Makefile`, GitHub Actions.
2. The five flat scripts moved from `python_scripts/` into an installable,
   typed package at `src/ms_denovo_db_utils/`, with byte-identical output.
3. Two order-dependence defects fixed (§9).
4. Full unit and container test suites.

### What has not been done

See §10 for the full list. The three that matter most:

- **CI has never been observed to pass.** The workflows were written and pushed
  but nothing has confirmed a green run. Check the Actions tab.
- **Nothing publishes.** `release.yml` builds and verifies only; quay.io is not
  linked to this repository.
- **The fixes have never been measured against real data.** They are proven
  correct on constructed fixtures and provably change nothing on well-behaved
  input, but nobody knows how many peptides they move on a real experiment.

---

## 3. The science, in brief

The pipeline searches the same spectra two ways and reconciles the results
against a large annotated protein library:

```
spectra ──┬─→ Comet (database search)   ─→ comet_peptides.txt   ─┐
          └─→ Casanovo (de novo)        ─→ casanovo_peptides.txt ─┤
                                                                  ├─→ combined_results.fasta
                                                                  │        │
                          annotated library + reversed decoys      │        ↓
                                        │                          │   DIAMOND blastp
                                        └──────────────────────────┴────→  │
                                                                           ↓
                                                                    reset_input.txt
                                                                           ↓
                                                                   percolator_RESET
                                                                     (FDR estimate)
```

### The important inversion

`build_reset_input` does **not** emit one row per PSM. A row is a **library
peptide region**: the subject subsequence that a DIAMOND alignment landed on.

- Every query peptide (from Comet or Casanovo) whose best alignment resolves to
  the *same* subject subsequence is aggregated onto *one* row.
- `Label` is `1` or `-1` according to whether the **annotated-library** protein
  carries the library decoy prefix — not according to Comet's own decoys.
- Consequently FDR is estimated over **annotated-database identifications**,
  not over spectrum identifications.

Two different query peptides landing on one library region is normal and is
exercised by the fixtures. Understanding this is the single most important
thing for working on `reset_input.py`.

### Two independent decoy universes

Easy to confuse; they are unrelated and have different defaults in the
pipeline:

| Prefix | Default | Applies to | Consumed by |
| --- | --- | --- | --- |
| `comet_decoy_prefix` | `COMET_DECOY_` | The small FASTA Comet searches | `process_comet_results` → `is_decoy` column |
| `library_decoy_prefix` | `LIBRARY_DECOY_` | The large annotated library | `build_reset_input` → `Label`, ambiguity filter |

A Comet decoy hit whose group resolves to a library **target** is contradictory
evidence and is discarded with a warning. A Comet decoy hit on a library
**decoy** is the consistent case and is kept.

---

## 4. Repository layout

```
├── src/ms_denovo_db_utils/     the package (all logic lives here)
│   ├── massutil.py             mass arithmetic, ppm error, rank scoring
│   ├── comet.py                Comet .txt  → peptide table
│   ├── casanovo.py             Casanovo mzTab → peptide table
│   ├── collate.py              peptide tables → query FASTA
│   ├── fasta.py                FASTA I/O (gzip-aware), reversed decoys
│   ├── diamond.py              DIAMOND outfmt 6, subject-region resolution
│   ├── reset_input.py          the percolator_RESET feature table
│   └── cli/                    one argparse entry point per tool
├── docker_bin/                 3-line shims → /usr/local/bin/*.py
├── tests/
│   ├── canonical.py            normalises environment noise before comparing
│   ├── unit/                   module behaviour, golden, determinism, perf
│   ├── docker/                 image contract + golden through the image
│   └── fixtures/               inputs; see tests/fixtures/README.md
│       └── golden/             expected outputs — see §8
├── Dockerfile, entrypoint.sh
├── Makefile                    the dev harness; CI runs the same targets
└── pyproject.toml              packaging + ruff/mypy/pytest/coverage config
```

`docker_bin/` exists solely to preserve the `/usr/local/bin/*.py` paths. Each
shim imports `main` from the matching `cli` module and calls it. The console
scripts declared in `pyproject.toml` (`process-comet-results`, etc.) expose the
same functions under clean names so the pipeline can migrate later.

---

## 5. The five tools

Argument order is the deployed contract. Do not reorder positionals.

### `process_comet_results.py [--decoy_prefix P] FILE...`

Collapses Comet PSMs to the best row per distinct peptide, lowest e-value wins.
Aggregates spectrum and peptidoform counts over **all** PSMs of that peptide,
across all input files.

Output columns (tab-delimited, to stdout):

```
plain_peptide  charge  e-value  protein  file  tryptic_n  tryptic_c
num_spectra  mz_ppm_error  is_decoy  proteins  rank_score  num_peptidoforms
```

- `tryptic_n` reads the residue *preceding* the peptide, which Comet embeds in
  its `modified_peptide` column as `K.PEPTIDER.A`. `-` (protein terminus)
  counts as tryptic.
- `tryptic_c` is `1` when the peptide ends in R/K **or** sits at a protein
  C-terminus (`modified_peptide` ends in `-`).
- `is_decoy` is `1` only when **every** protein in the comma-separated list
  carries the prefix. A peptide shared with any target is not a decoy.
- `proteins` duplicates `protein`. This is deliberate; see §7.

### `process_casanovo_results.py FILE...`

The same, for Casanovo mzTab, keyed on the **highest** `search_engine_score[1]`.

```
peptide_sequence  charge  search_engine_score[1]  file
mz_ppm_error  num_spectra  rank_score  num_peptidoforms
```

- The sequence has inline modification masses stripped (`[^A-Z]` removed), so
  `M+15.995DLGEEHFK` → `MDLGEEHFK`. The unstripped form is the peptidoform.
- **Negative scores get +1.** Casanovo subtracts 1 when a PSM fails its
  precursor m/z filter; adding it back deliberately disables that filter,
  because the homology search — not Casanovo — decides what is credible.

### `collate_into_fasta.py FILE...`

Union of column 1 of each file (header skipped) as a FASTA where each peptide
is its own entry named after itself. DIAMOND carries the query name through to
its output, which is what lets results be joined back. Output is **sorted**.

### `generate_reverse_decoys.py INPUT [--decoy_prefix P]`

Writes every target entry followed by its reversed decoy. Gzip is detected by
**magic number, not extension**; a gzipped input produces gzipped stdout, and
`GENERATE_LIBRARY_DECOYS` in the pipeline names its output `.fasta.gz` or
`.fasta` on exactly that behaviour. A trailing `*` is stripped before reversal.

### `build_reset_input.py COMET CASANOVO DIAMOND FASTA LIB_PREFIX [COMET_PREFIX]`

The feature table for percolator_RESET. `FASTA` is the annotated library
**with** decoys. `COMET_PREFIX` is accepted and ignored (§7).

21 columns:

```
SpecId  Label  ScanNr  database_peptide_length  max_diamond_bitscore
max_diamond_perc_identity  num_casanovo_peptides  num_comet_peptides
casanovo_num_spectra  casanovo_best_score  casanovo_ppm_error
casanovo_num_peptidoforms  comet_num_spectra  comet_n_tryptic  comet_c_tryptic
comet_best_score  comet_ppm_error  comet_num_peptidoforms
combined_rank_score  Peptide  Proteins
```

- `SpecId` is `<library peptide>_1`; `ScanNr` counts from 1 in output order.
- `comet_best_score` is `log10(1 + 1/(e-value + 1e-20))` — larger is better,
  and the floor exists because Comet does report an e-value of exactly `0`.
- `combined_rank_score` is `4 − comet_rank − casanovo_rank`; an engine that
  contributed nothing supplies the worst rank, `2.0`.
- Rows where every Comet hit was discarded and Casanovo contributed nothing are
  omitted entirely rather than emitted with zero evidence.
- Peptides absent from the DIAMOND results are skipped with a sorted warning on
  stderr.

---

## 6. Invariants

These are enforced by tests. Breaking one is a behaviour change, not a
refactor.

1. **Output is deterministic.** Identical inputs produce byte-identical output
   regardless of `PYTHONHASHSEED`. Any new code that iterates a set, or a dict
   built from a set, and affects output must impose an order. `tests/unit/
   test_determinism.py` runs the CLIs across eight seeds. **Sort by default.**
2. **The `/usr/local/bin/*.py` paths and their argument orders are fixed.**
3. **A gzipped library in means a gzipped library out.**
4. **The annotated library is streamed, never loaded whole.** It can be very
   large. `add_subject_sequences` reads one protein at a time on purpose.
5. **Malformed input is an error with a non-zero exit.** Nextflow reads exit
   codes; the original scripts printed a message and exited 0, producing an
   empty result file that read as success.
6. **The image keeps `bash` and `procps`.** Every nf module writes through
   `> >(tee ...)` process substitution, which needs real bash; `ps` backs
   Nextflow's task metrics.
7. **The entrypoint passes commands through and propagates exit codes.**

---

## 7. Deliberate decisions that look like bugs

Do not "fix" these without reading the reasoning.

**The duplicated `proteins` column.** `comet_peptides.txt` emits the protein
list twice, as `protein` and `proteins`. Nothing in `nf-ms-denovo-db` reads it,
but the file is published to `results/` and consumers outside these two
repositories cannot be ruled out. The gain from removing it is cosmetic; the
risk is a silent break in somebody's downstream tooling.

**`build_reset_input` accepts an ignored `comet_decoy_prefix`.** Comet decoy
status is already resolved upstream into the `is_decoy` column. The positional
remains optional-and-ignored so a pipeline pinned to an older revision keeps
working. Removing it requires a coordinated change in both repos.

**Some fields are carried as `str`, not parsed.** `mz_ppm_error`, `tryptic_n`
and `tryptic_c` pass from the intermediate tables straight into
`reset_input.txt` as text, and several numeric accumulators are initialised to
int `0` rather than `0.0`. This is to preserve output formatting exactly:
`str(0)` is `"0"` but `str(0.0)` is `"0.0"`, and `0.40` re-parsed and
re-rendered becomes `0.4`. Changing these types changes the bytes of a
published result file. The golden tests will catch it — understand it before
overriding it.

**`write_reset_input` resolves `sys.stdout` inside the function.** Not as a
default argument, which would bind whatever the stream was at import time.

---

## 8. Development

Requires Python 3.10+ and, for container tests, Docker.

```sh
make install      # create .venv, install the package and dev tools
make check        # everything CI runs: lint, typecheck, test, image, container tests
```

| Target | Does |
| --- | --- |
| `lint` | `ruff check`, `ruff format --check`, `hadolint` if installed |
| `format` | `ruff check --fix` and `ruff format` |
| `typecheck` | `mypy --strict` over `src/` and `tests/` |
| `test` | unit tests with coverage; floor 95% |
| `image` | `docker build -t ms-denovo-db-utils:dev .` |
| `test-docker` | rebuilds the image, then runs container tests |
| `clean` | remove `.venv` and caches |

CI (`.github/workflows/ci.yml`) runs four jobs: lint, typecheck, a test matrix
over 3.10–3.13, and an image build plus container tests. **3.10 is the
contract** — it is what `python:3.10-slim` ships — and the newer legs are
`continue-on-error` advisory warnings for a future base-image bump.

### Test layout

| Path | Covers |
| --- | --- |
| `tests/unit/test_massutil.py` | ppm error incl. 13C tolerance, rank scoring |
| `tests/unit/test_comet.py` | trypticity, decoy status, PSM collapsing |
| `tests/unit/test_casanovo.py` | mod stripping, score lift, collapsing |
| `tests/unit/test_fasta.py` | gzip, parsing, decoy generation |
| `tests/unit/test_diamond.py` | parsing, tie-breaks, ambiguity, subject regions |
| `tests/unit/test_reset_input.py` | grouping, labels, aggregation, ordering |
| `tests/unit/test_golden.py` | full chain vs. stored expected output |
| `tests/unit/test_determinism.py` | eight hash seeds, in subprocesses |
| `tests/unit/test_performance.py` | rank scoring on 200k peptides (`slow`) |
| `tests/unit/test_packaging.py` | version metadata, release gate inputs |
| `tests/docker/test_image_contract.py` | script paths, entrypoint, bash, `ps` |
| `tests/docker/test_pipeline_output.py` | golden output *through the image* |

### How the golden files work

`tests/fixtures/golden/` holds output captured from the **original** flat
scripts before any restructuring. That is what proves the move into `src/`
preserved behaviour, and what pins everything since.

> **A diff in a golden file is a change in scientific output.** Never
> regenerate one to make a test pass. Work out why it changed, decide whether
> the change is correct, and if so update the golden in the same commit that
> causes it, with the reasoning in the commit message.

The same goldens are asserted by the unit tests and, through the built image,
by the container tests. `tests/canonical.py` normalises only the input-path
column, which legitimately differs between a local run and a container run.

---

## 9. History: the two defects that were fixed

Useful context, because the fixture set is shaped around them.

**Order-dependent results.** `build_reset_input` computed `best_diamond_label`
while iterating a set and read that same label to decide whether to discard a
decoy Comet hit. Python randomises string hashing per process, so the outcome
varied run to run. Measured on the committed fixtures: three of eight seeds
aborted with a `ValueError` and five succeeded; on another fixture seed 0
counted a decoy Comet hit into a target row (`comet_num_spectra=2`) while seed
1 discarded it (`=1`). Row order and `ScanNr` varied on every dataset.

Fixed by resolving the group's best hit *before* examining any Comet hit, and
ordering groups and members. The `ValueError` guard was removed: it fired only
because the label was still being computed while those hits were consumed, and
it aborted whole pipeline runs.

**Order-dependent ambiguity filter.** `read_hits` cleared its accumulated
protein set whenever a better e-value arrived, so a peptide aligning to both a
target and a decoy was only filtered when the rows happened to arrive in one
order. Fixed by accumulating across every alignment before checking. Best-hit
tie-breaking is now by descending bit score then subject name, rather than by
arrival order.

**Result:** on well-formed data nothing changed. `reset_input.txt` for the main
fixture dataset is byte-identical to the pre-fix golden and now identical
across every seed tested. Only the pathological cases behave differently.

---

## 10. TODOs and open questions

**Verify CI actually passes.** The workflows were pushed but no run has been
observed. Locally verified as proxies: the full suite passes inside
`python:3.10-slim` (Python 3.10.20, including the `tomli` fallback path that
the 3.12 leg skips), `hadolint` exits 0, the image builds, and the release
version-gate logic behaves. The YAML itself is unproven; expect a first-run
fixup.

**Measure the fixes against real data.** The highest-value outstanding task.
Run `1.1.10` and this revision over the same real experiment and diff
`reset_input.txt`. The fixtures prove the fixes *correct*; only real data shows
*how much they move*. Nobody has this dataset yet.

**Wire up publishing.** `release.yml` is build-and-verify only. When quay.io is
linked, add `QUAY_USERNAME` / `QUAY_ROBOT_TOKEN` secrets and the login/push
steps — they are written out in a comment at the bottom of that file.

**Release `1.2.0`.** The version is `1.2.0.dev0`. Bump it in `pyproject.toml`
before tagging, or the release workflow will reject the tag: it refuses to
proceed when the git tag and the declared version disagree. That check, plus
the `org.opencontainers.image.revision` label stamped from `VCS_REF`, is what
makes a running image traceable to a commit — previously impossible, since the
tag existed only as a string in the pipeline's config.

**Bump the pipeline.** `nf-ms-denovo-db/container_images.config` still pins
`1.1.10`. It must move once an image is published.

**Optional, needs a coordinated two-repo change:** drop the ignored
`comet_decoy_prefix` positional from both the tool and
`modules/build_reset_input.nf`.

**Optional, needs a decision:** remove the duplicated `proteins` column (§7).

---

## 11. Footguns

**A stale image gives false green.** `pytest -m docker` reuses an existing
`ms-denovo-db-utils:dev` and only builds when the tag is absent. After editing
`src/`, run `make test-docker` (which rebuilds) or `make image` first —
otherwise you are testing the previous build. Override the tag with
`MS_DENOVO_DB_UTILS_IMAGE`.

**`make check` on a machine without Docker is not full coverage.**
`test-docker` prints `SKIP` and exits 0 by design, so `make check` stays usable.
A green run there says nothing about the image.

**Hash randomisation is the recurring hazard.** See invariant 1. Sort anything
that reaches output.

**`tomllib` is 3.11+.** The 3.10 leg — the one that matters, since it is the
image interpreter — uses the `tomli` backport, declared as an environment-marked
dev dependency. Code touching TOML needs the version-guarded import already in
`tests/unit/test_packaging.py`.

**Tests are type-checked too.** `mypy --strict` covers `tests/`; new test
functions need `-> None`. An unannotated test can silently stop asserting what
it claims to.

**The pipeline repo can be pinned to an old image.** Any change to a CLI
contract must keep working for a pipeline revision that has not been updated.
Additive and optional, never reordered or removed.
