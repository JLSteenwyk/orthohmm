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
packaging tests. Biopython, Pillow and setuptools have subsequent test-only
security updates; other pins retain the recorded local choices. They are not
substituted into historical scientific environments. A clean
test environment needs these packages even when optional application features
or publication workflows are not being executed.
NumPy is pinned to 2.2.6 because retained factorial, sequence-control and
matched-graph numerical replay compares the engine version as well as the
statistics. This preserves exact provenance comparisons; it does not loosen
their equality checks or change the application's broader runtime requirement.
The [security follow-up](../benchmark_tools/results/TEST_DEPENDENCY_SECURITY_20261001.md)
retains fifteen test-manifest advisory ranges, regression checks against old
and current pins, and fresh private-environment validation. Range exclusion
does not establish complete security or clear historical benchmark locks.

Actual affinity/cgroup observations require Linux facilities. Importing the
observer and testing it with injected readers does not require those facilities;
requesting native observation without the affinity API fails explicitly.
Native calibration fixtures report capability skips on unsupported platforms.
Synthetic release-guard and environmental-review cases inject a fixed test-only
boot identity rather than reading the host's Linux boot file. Other file reads
and real native identity tests remain unchanged. Short assessment fixtures
resolve the temporary root (including macOS's `/tmp` symlink); the production
admission gate still rejects indirect destinations and mismatched boot identities.
Synthetic measurement, evidence-replay, cgroup-frontier and service-reply tests
also explicitly request that fixture. Missing boot data and reboot rejection
remain tested; production collectors still require real Linux capabilities.
See the [scoped validation](../benchmark_tools/results/CI_SYNTHETIC_CLOCK_FIXTURES_20261001.md).
Historical QfO handoff tests explicitly relocate only root/interpreter bindings
in fresh temporary copies, exercising dummy executors and synthetic scheduler
rows with shell-special paths. Archived batches and their provenance checks
stay unchanged. See [validation and production boundaries](../benchmark_tools/results/CI_QFO_BATCH_FIXTURES_20261001.md);
test portability does not make historical native workflows portable.
Batch handoff checks require a Bash runtime that stops after failed standalone
guards under `set -e`; direct semantic tests do not skip unsupported behavior.
The macOS test jobs prepare an explicit Homebrew Bash on PATH and log its version
alongside the system version. The [observed failure and amendment](../benchmark_tools/results/CI_BASH_GUARD_RUNTIME_20261001.md)
distinguish the new setup from unverified remote execution and historical runtime
claims. Archived scripts must not be assumed safe under an arbitrary shell.
Synthetic affinity cases install their injected APIs explicitly on unsupported
hosts; genuine native workload/resource/library-worker cases retain capability
skips rather than mocked native results. Fake-service capture hashes its own
temporary declared file. The [capability validation](../benchmark_tools/results/CI_AFFINITY_FIXTURES_20261001.md)
distinguishes actual Linux execution from a deliberate API-absence simulation.

Pipeline metrics export is tested with and without Linux `/proc`: absent RSS
retains the explicit `unavailable` convention and zero sentinels, not measured
zero memory. Only the positive native RSS observation requires `/proc`;
synthetic readers, export and rejection cases still run. Both OrthoBench resource
readers reject unavailable memory. See the
[capability checks and sampling limits](../benchmark_tools/results/CI_RSS_CAPABILITY_TESTS_20261001.md).
The [subsequent inspected macOS logs](../benchmark_tools/results/CI_NATIVE_CHECKOUT_SETUP_20261001.md)
confirm the scoped Bash/affinity/RSS corrections, not broader CI success.
MacOS test jobs now prepare a job-local GNU GCC driver and explicitly build
all three baseline CPU libraries for checkout imports as well as installation.
The checkout helper rejects inherited, incomplete, unloadable or nonbaseline
libraries before exposing them. This prospective native setup still requires
remote verification; it does not change scientific defaults or benchmark builds.
Synthetic long-worker tests explicitly inject cgroup identity and resource
snapshots; their short subprocess checks do not poll host counters. Actual
`/dev/shm` copy cases skip when unavailable, without substituting disk storage.
Portable path guards still test all three layouts and reject unsafe input/output
bindings and pre-phylogeny commands. The
[scope and tmpfs validation](../benchmark_tools/results/CI_CGROUP_TMPFS_FIXTURES_20261001.md)
distinguishes fixture evidence from real preparation or memory charging.

Unit, fast and unit-coverage targets include top-level `tests/test_*.py` as
well as `tests/unit`; integration remains a separate target. Fast tests exclude
the `slow` marker, not the top-level regression cases. CI keeps every Python
matrix job's outcome rather than cancelling siblings after the first failure.
The two test jobs fetch full Git history because frozen-source archive tests
read an actual historical scientific commit, not HEAD. Matched-graph rendering
still requires exact replay and now reports differing JSON paths/values on
rejection; one-ULP and engine-metadata changes are not silently tolerated.
An explicit `count-level` rendering policy reuses the pre-existing independent
1e-12 count-reproduction contract and records all accepted differences; it
does not replace the default strict guard or relax metadata. The retained
plot test chooses that portable mode; strict behavior has a separate
current-runtime fixture and one-ULP rejection tests. See the
[policy and evidence](../benchmark_tools/results/PORTABLE_MATCHED_GRAPH_RENDERING_20261001.md).
The [subsequent macOS log inspection](../benchmark_tools/results/CI_SYNTHETIC_CLOCK_FIXTURES_20261001.md)
confirms those 61 rendering and seven count-reproduction cases pass. The broader
CI suite still fails; this is not native scientific or full-suite confirmation.

A separate Linux CPU-wheel job builds committed package sources in a clean
directory, installs into a fresh private venv and runs the existing isolated
installed-package verifier. It exercises builtin/Leiden standard and
high-sensitivity fixtures, not the full phylogenetic pipeline or benchmark
accuracy. Existing macOS test jobs are retained. Result and command logs are
uploaded even after normal failures; this is not a replacement for the full
suite or proof of a hermetic runtime. The [execution report](../benchmark_tools/results/CI_CPU_WHEEL_INSTALLATION_20261001.md)
distinguishes the local installation attempt. The [inspected remote result](../benchmark_tools/results/CI_REMOTE_CPU_WHEEL_RESULT_20261001.md)
confirms the separate Linux job passed; it does not clear the macOS matrix.

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
