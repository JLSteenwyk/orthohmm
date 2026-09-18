# Tests

Run unit tests with`make test.unit` and native integration tests with
`make test.integration`. Pytest discovery is limited to`tests/` so retained
benchmark worktrees are not collected. Integration outputs and copied input
FASTAs live under pytest temporary directories; the commands do not remove
or overwrite`tests/samples` output files.

Integration tests require working`phmmer` and`mcl` executables. The MCL build
must support the application's ABC-input and threading options. A legacy
OrthoMCL-specific installation should not replace the integration dependency.
Select executables for tests without changing the global PATH:

```sh
ORTHOHMM_TEST_MCL=/absolute/path/to/mcl \
ORTHOHMM_TEST_PHMMER=/absolute/path/to/phmmer \
make test.integration
```

The simple and long-name fixtures run with both forward and reversed input
creation orders. Checks cover complete reference group memberships, every
species count, every output FASTA sequence, single-copy occupancy, exact
file inventories and species-qualified single-copy headers. Group numbering,
species-column order and FASTA record order are not biological invariants.
Deliberate corruption tests verify that semantic errors still fail.

Historical single-copy text fixtures list all groups and are not a valid
oracle for single-copy selection. The current test independently derives
that selection from the unchanged expected memberships, input species and
the configured occupancy threshold; the emitted ID list must agree with the
single-copy FASTA directory. The historical files remain unchanged to retain
the evidence of the old output bug.
