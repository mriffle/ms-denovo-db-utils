# Test fixtures

Small, hand-checkable inputs in the real formats the pipeline consumes. The
library proteins are built by concatenating known peptides, so every DIAMOND
`sstart`/`send` can be verified by inspection.

## Main dataset

`comet/run{1,2}.txt`, `casanovo/run{1,2}.mztab`, `diamond/hits.dmnd.txt` and
`fasta/library_targets.fasta` form one coherent mock experiment covering:

- a peptide found by both engines, across two runs, with several peptidoforms
- Casanovo-only and Comet-only peptides
- a precursor misassigned by one 13C, which the ppm logic must recover
- tied Comet e-values (rank ties) and an e-value of exactly zero
- protein N-terminal (`-.PEP.X`) and C-terminal (`X.PEP.-`) peptides
- a Comet decoy hit, and a peptide mapping to both target and decoy proteins
- two distinct query peptides landing on the same library subsequence
- peptides absent from the DIAMOND results (missing-peptide warning)

## Pathological datasets

Isolated so the main dataset stays deterministic:

- `diamond/ambiguous.dmnd.txt` — one peptide with its target hit listed first
  and one with its decoy hit listed first. Both are equally ambiguous and must
  be filtered identically.
- `comet/crash.txt` + `diamond/crash.dmnd.txt` + `fasta/crash_plusdecoys.fasta`
  — a group containing both a library-decoy and a library-target hit, where a
  decoy Comet hit is the group's best. The library decoy region and the target
  protein share a subsequence, so both queries land in one group.
- `comet/decoy_meets_target.txt` + `diamond/decoy_meets_target.dmnd.txt` — a
  decoy Comet hit landing on a library target subsequence.

`fasta/with_stop_codon.fasta` pins the trailing `*` stripping in decoy
generation.
