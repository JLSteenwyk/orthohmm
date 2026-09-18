# Single-Copy List And Integration Repair

The investigation of[retained failures](REGRESSION_AUDIT_970b496_20260918.md)
identified both order-sensitive tests and a real output bug. The single-copy
ID writer iterated every copy-number-table key rather than the selected
single-copy group list. Its text file therefore disagreed with the single-copy
FASTA directory. The caller now passes the same selected IDs used by that
directory writer; an empty selection produces an empty list file. IDs match
the FASTA stems, including the existing unpadded stem convention.

This changes an output summary only. No search, grouping, reconciliation,
threshold, or biological selection algorithm was changed. Existing frozen
benchmark executors remain at their original commits; their orthology
predictions and endpoint scores are not regenerated or silently replaced.

## Stronger Integration Checks

Replaced redundant native invocations and byte-order assertions with two
fixtures each run in forward/reversed input creation order. Every run checks:

- The full unchanged expected partition:38simple or1155long-name genes.
- Every copy-number row by species and actual group membership.
- Every group FASTA identifier and full sequence against original inputs.
- Complete output file inventories:5/988all-group files and1/9single-copy files.
- Independent single-copy selection from species uniqueness and strict
  occupancy greater than0.5, including the exact ID list and prefixed headers.

All files are generated in pytest temporary directories. Expected biological
partitions and historical output fixtures were not rewritten. The old
single-copy text fixtures are explicitly documented as invalid selection
oracles. Corruption tests reject altered counts, duplicate species columns,
changed sequences, duplicate FASTA records, wrong species prefixes, extra or
missing single-copy IDs, extra output files and changed/duplicated membership.

Pytest discovery now targets`tests/`, not retained benchmark worktrees.
Makefile commands explicitly select unit/integration directories and no
longer delete sample outputs. Scoped executable overrides are documented in
`tests/README.md`; this host requiresMCL14-137 rather than the shell's legacy
OrthoMCL-specific02-063 binary for these integration tests.

## Verification

-60focused unit/helper/corruption tests pass.
-Full unit suite:5272passed,9skipped in91.91seconds.
-Four native integration cases pass in70.41seconds withMCL14-137 and the
  existing phmmer dependency. This replaces eight mostly redundant legacy
  cases with four complete-output runs, not four selected assertions.
-Makefile dry-run and scoped whitespace checks pass.

[Integration JUnit](integration_isolated_outputs_20260918.xml) SHA-256:
`6fcb93fc73d7ebdfc63664a32c8f631fa5a822230266c67acd39ec9c7eea8c4d`.
Full unit JUnit retained at`benchmarks/work/unit_output_fix_20260918.xml`,
SHA-256`3738d391ef5ba172e72d55f82507863e3c6cdd2099e7f93a0c92d94709c81fc2`.
These are regression results, not evidence of benchmark superiority or
publication readiness. Prior generated sample changes remain unstaged.
