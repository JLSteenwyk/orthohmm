# Explicit Native Utilities For CI Tests

The inspected macOS unit run confirms the preceding native-build and
cgroup/tmpfs corrections, but still fails. GNU Time and Pandoc failures motivate
explicit test-job utility installation and narrowly scoped smoke-test binding.
Frozen production commands, source pins, scientific settings and resource
admission rules remain unchanged. The [machine-readable receipt](ci_native_utilities_20261002.json)
pins source files, native binaries, both local test reports and the remote log.

## Observed Remote Evidence

Source `42e48844` Python 3.11 fast job `110690397294` reports 13,665 passed,
60 failed, four errors, 112 skipped and 30 warnings in 544.21 seconds. Its log
records actual Homebrew GCC 16.2.0, test Bash 5.3.15 and system Bash 3.2.57;
the fresh checkout builder reports three load-verified baseline CPU libraries.
Within the previous 70-case native scope, 67 pass and three skip. Scalar/Numba
parity, profile expansion, fresh staged compilation and child alignment pass;
AVX2, CUDA and the separate baseline wheel compile case skip. The last skip's
specific cause has not been isolated. These are bounded test observations,
not complete native/platform compatibility or scientific benchmarking.

Within the previous 59-case cgroup/tmpfs scope, 51 pass and eight actual
tmpfs-copy cases skip. Synthetic orchestration, unavailable-cgroup rejection
and portable layout guards execute. No disk-backed tmpfs substitute or resource
admission is added.

The same log shows BSD `/usr/bin/time` rejecting GNU options, absent timing
logs in runner smoke tests, and the real citation parser failing because Pandoc
is missing. At 03:39:06 UTC on 2 October all five source-42e test jobs are
terminal failure; wheel/docs succeed. Only the Python 3.11 fast log is inspected
here. Do not infer sibling failure causes, full integration execution or success.

## Scoped Utility Binding

Only the two macOS test definitions gain one native-utilities setup step. It
installs [Homebrew GNU Time](https://formulae.brew.sh/formula/gnu-time) and
[Pandoc using its documented macOS installation](https://pandoc.org/installing.html),
records actual versions, and exports the absolute `gtime` path as
`ORTHOHMM_TEST_GNU_TIME`. Existing steps and other jobs are structurally unchanged.
No system binary is replaced, and no utility is installed on the shared host.
The Homebrew formula's advertised version is not a verified future CI version.

The explicit GNU Time fixture validates an absolute executable and successful
GNU version output. Invalid configured paths, BSD binaries, probe errors and
timeouts fail closed without searching for a replacement. When unconfigured,
discovery checks `gtime`, `time` and `/usr/bin/time`; missing GNU capability
skips only explicitly requested native smoke cases. Parser, path, allocation,
early provenance and other independent guard cases still run.

The companion test uses its existing explicit executable argument. Seven runner
modules are bound only inside requested tests: a module-local subprocess proxy
replaces a leading literal `/usr/bin/time` executor, preserving all subsequent
arguments, environments, output handling and native child execution. It does not
mutate original argv or patch the global subprocess module. Teardown restores
the module namespace. Historical command metadata is not reinterpreted: this
binding is for temporary smoke tests, not frozen production execution.

## Executed Validation

**100 focused Linux cases pass in 3.10 seconds**, zero failures/errors/skips,
including 24 new selection/binding tests. Real GNU Time runs retain child
exit 0/7, native output and verbose timing logs. Existing tests cover owned
timeout cleanup, retained failures, stop/continue rules, postflight provenance,
no implicit restart and explicit pending-admission states. The real Pandoc
parser test also passes using the existing local Pandoc 3.1.3. Local GNU Time
reports version `UNKNOWN`; its binary and complete version output are pinned
rather than assigned an assumed release number.

A separate process deliberately makes GNU version probes fail as BSD Time
would. The same 100-case scope gives **76 passes and 24 explicit native smoke
skips in 2.29 seconds**, zero failures/errors. This is a local capability
simulation, not actual macOS execution; the two panels are not additive.
Explicit invalid configuration tests still execute and reject.

Structured PyYAML 6.0.1 comparison proves exactly one added utilities step per
test job, with all existing steps and other jobs unchanged. All eight production
utility/runner files are byte-equal to the base commit. Both JUnit files and
twelve changed workflow/test source pins are retained in the receipt.

## Remaining Work

Commit/push and inspect the new automatic CI run; utility installation and
new native smoke outcomes on macOS remain unverified. Workstation-only paths,
missing retained raw data, other GNU utility/Linux runtime dependencies and
remaining test failures still need separate assessment. Tiny command timings
are not controlled scientific resources; GNU Time's maximum process RSS is
not a simultaneous process-tree or cgroup peak.

No scientific input/default/score, production fallback or admission changes,
shared-package upgrades, local workload polls/questions, DGX work or unrelated
job/service actions occur. Controlled timing stays deferred. Remaining QfO
uncertainty, source/rights review, comparable resources and the final executable
versioned release/archive prevent publication readiness; the full goal stays
active.
