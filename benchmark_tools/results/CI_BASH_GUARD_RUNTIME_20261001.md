# Bash Runtime Guard Validation

The latest inspected remote handoff tests exposed fail-open guard behavior,
not a successful portability correction. Require direct shell semantic checks
and prepare an explicitly installed test Bash runtime without changing archived
scientific batch scripts. Continue from pushed `9674fe93`.

## Observed Failure

Source-bec79e1d macOS Python 3.12 job 110674083803 reports **13,458 passed,
198 failed, 28 errors, 96 skipped, 30 warnings in 411.44 seconds**. All eight
new relocation-fixture checks and all 41 reconciliation-admission cases pass.
Other handoff cases now expose success returns for invalid commits, commit/job
formats, parity and task indices. Prior local checks did not demonstrate remote
guard correctness. Preserve these failures and their
[log receipt](ci_bash_guard_runtime_20261001.json); do not change their expected
rejection or classify them as passing. The actual failing Bash version was not
recorded, so the shell-version explanation remains a hypothesis, not a diagnosis.

At 02:31:42 UTC source-bec run 36954409877 is terminal failure across all five
test jobs, with CPU-wheel/docs success. Only Python 3.12's log is inspected.
Source-9674 run 36955273317 has CPU-wheel/docs success and five live test jobs
at the subsequent checkpoint; do not infer their outcomes or restart either run.

## Prospective Runtime And Checks

Add six direct Bash tests: valid-guard continuation and failure stops for a
simple command, commit equality, commit format, job format and task parity.
They require nonzero failure with no downstream marker, and report selected
shell path/version on rejection. No skips, expected-return relaxation or changes
to the archived guard expressions are introduced. The semantic contract follows
the [GNU Bash set documentation](https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html).

Only the two macOS test-job definitions gain a Bash installation step using
the [Homebrew Bash formula](https://formulae.brew.sh/formula/bash). Record both
its version and `/bin/bash`'s version, and prepend its bin directory using
[GitHub's job PATH mechanism](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-commands#adding-a-system-path).
The next automatic CI must verify the actual selected shell and guard outcomes.
This is a prospective test-runtime amendment, not proof of its execution or a
retroactive alteration of historical environments. Homebrew's resolved package
is not a frozen scientific-runtime pin or comprehensive compatibility guarantee.

**261 focused cases pass in 12.89 seconds**, zero failures/errors/skips, including
six new guard cases and all existing handoff/replay/relocation checks. Actual
local Bash is 5.2.21; its bytes and JUnit are pinned in the receipt. Structured
YAML comparison proves exactly two added test setup steps, with all other jobs
and settings unchanged. That read-only check used existing PyYAML 6.0.1 after
the private test environment lacked YAML; no package was installed or upgraded.
All eight archived batch scripts remain byte-equal to the prior commit and their
existing receipt pins validate. Local and remote panels are not additive.

Commit/push, then inspect actual new CI without retries. RSS capability tests
remain a separate follow-up: existing collectors already label unavailable
Linux `/proc` memory, and benchmark readers reject that convention. No metrics
code or tests are edited in this milestone. Historical guard behavior is not
claimed portable to an arbitrary Bash runtime; future release execution must
bind and validate its shell rather than relying on a generic executable name.
Scientific inputs/defaults/scores, native engines and historical results remain
unchanged. No new native benchmark, controlled timing, host poll/question, DGX,
shared-package or unrelated workload/service action occurs. Broader CI,
source/rights, uncertainty, comparable resources and public release remain open.
