# CI Execution And Fixture Portability

The [CI run at ea90d57f](https://github.com/JLSteenwyk/orthohmm/actions/runs/36941807688)
is terminal failure. Dependency installation and collection now succeed;
documentation succeeds, but all five test jobs fail during execution. This is
not a passing clean installation or a reclassification of the earlier failures.
The [measured receipt](ci_fixture_portability_20261001.json) pins terminal API
observations, all five downloaded test logs, current fixtures and local JUnit.
Credentials and signed redirect URLs were not retained or forwarded to log hosts.

| Job | Passed | Failed | Errors | Skipped | Test duration (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fast Python 3.10 | 13,208 | 265 | 109 | 94 | 616.22 |
| Fast Python 3.11 | 13,204 | 269 | 109 | 94 | 645.00 |
| Fast Python 3.12 | 13,201 | 272 | 109 | 94 | 531.06 |
| Fast Python 3.13 | 13,193 | 280 | 109 | 94 | 463.50 |
| Full Python 3.11, unit coverage stage | 13,204 | 269 | 109 | 94 | 716.67 |

Each invocation accounts for 13,676 unit/top-level cases. Failed Make stages
prevented their subsequent integration stage; the durations are CI test times,
not scientific benchmark runtimes or controlled timing evidence.

## Focused Corrections

The short private-assessment fixture now resolves its temporary root. macOS
resolves `/tmp` to `/private/tmp`; the old unresolved path was rightly rejected
by the production admission gate. A new symlink-root regression exercises this
on Linux too. The gate's direct-path, symlink, no-overwrite and allocation checks
remain unchanged; existing rejection tests are retained.

Synthetic release-guard/environmental-review fixtures now inject a fixed
test-only boot identity rather than requiring the host's `/proc` boot file.
The narrowly scoped fixture preserves every other file read. Added checks
cover the fixture boundary and a boot change during review: review fails and
does not create a release signal. Real loaded-executable identity checks remain
enabled on Linux and explicitly skip on non-Linux systems; synthetic workflow
tests are not skipped. No live host observation or timing admission is fabricated.

The three production admission/execution/worker files are byte-equal to HEAD.
No benchmark score, frozen method, scientific default, source binding or
historical receipt is rewritten by these fixture changes.

## Executed Validation

In the previously prepared private Python 3.12 dependency environment,
**225 focused cases pass in 6.73 seconds**, with zero failures/errors/skips.
This includes the fixture corrections, native Linux identity checks, observer
capability tests and dependency/discovery checks. Four new cases are included.
An earlier 158-case pass preceded the final native capability marker and is
not counted as 158 additional unique tests. The final JUnit binds the current
source snapshot. Scoped whitespace checks pass.

This is a focused Linux checkout result using retained native libraries, not
a full new regression, macOS confirmation, clean installed-package restoration
or successful remote CI. Observe the next pushed source's actual CI handle;
do not resubmit the failed run or infer remote success from these local checks.

## Remaining Work

Failures also include workstation-specific paths, absent external benchmark
assets, genuinely Linux-only affinity/cgroup/tmpfs tests, missing native kernels
and Apple clang's rejected OpenMP option. Numerical/provenance reproduction
failures remain to be traced individually. Do not blanket-skip test modules,
weaken safety gates, invent native capabilities, globally rewrite historical
paths or replace frozen benchmark dependencies merely to make CI green.

Controlled Threadripper timing stays deferred. No quiet-window question or
contention poll, DGX access, scientific rerun, shared-package upgrade or
unrelated job/service action occurred. The publication goal remains incomplete.
