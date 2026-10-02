# Tests

After installing OrthoHMM and its dependencies in a development environment,
install the test-only dependencies with:

```sh
python -m pip install -r tests/requirements.txt
```

This includes the pinned SQL parser used to inspect historical TreeFam and
Selectome archives without executing their SQL. Without `sqlglot`, two unit
modules are skipped during collection, hiding 13 test cases from the run.
Both CI test jobs install this file and check that the parser imports before
running tests. These are test dependencies, not application runtime
requirements or a complete lock for the publication benchmark environments.
The file also pins libraries used by workflow, biology, plotting, PDF and
packaging tests. Their versions match the recorded 1 October local regression;
they are not substituted into historical scientific environments. A clean
test environment needs these packages even when optional application features
or publication workflows are not being executed.

Actual affinity/cgroup observations require Linux facilities. Importing the
observer and testing it with injected readers does not require those facilities;
requesting native observation without the affinity API fails explicitly.
Native calibration fixtures report capability skips on unsupported platforms.
Synthetic release-guard and environmental-review cases inject a fixed test-only
boot identity rather than reading the host's Linux boot file. Other file reads
and real native identity tests remain unchanged. Short assessment fixtures
resolve the temporary root (including macOS's `/tmp` symlink); the production
admission gate still rejects indirect destinations and mismatched boot identities.

Unit, fast and unit-coverage targets include top-level `tests/test_*.py` as
well as `tests/unit`; integration remains a separate target. Fast tests exclude
the `slow` marker, not the top-level regression cases. CI keeps every Python
matrix job's outcome rather than cancelling siblings after the first failure.

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

For a complete default-discovery regression run, use the same explicit
dependency selection and retain a machine-readable report:

```sh
ORTHOHMM_TEST_MCL=/absolute/path/to/compatible/mcl \
python -m pytest -q --junitxml=/absolute/path/to/regression.xml
```

Ten installed-runtime OrthoMCL checks are opt-in through
`ORTHOHMM_LEGACY_BLAST_SMOKE=1`. They require the retained native installations
and, for staged inference, the dedicated Python environment under
`benchmarks/work/orthomcl_python_env_20260918`. They use temporary fixtures,
not the full benchmark datasets. Enable them only where those prerequisites
are available; report skips and separately executed native checks explicitly.
The compatible MCL selected for OrthoHMM integration tests does not replace
the legacy MCL used by the frozen OrthoMCL benchmark workflow.

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
