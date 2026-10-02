# RSS Capability Test Validation

Memory reporting must distinguish unavailable Linux observations from measured
resource evidence. This milestone tests the existing contract without changing
the collector, benchmark admission, scientific defaults or retained scores.
Continue from pushed `8b48b8bf`; controlled timing remains deferred and no quiet
window is needed for these tests.

## Observed Failure And Correction

The retained source-bec79e1d macOS Python 3.12 log shows both original metrics
tests failing on `0 > 0`. Its overall result is 13,458 passed, 198 failed,
28 errors, 96 skipped and 30 warnings in 411.44 seconds. The
[receipt](ci_rss_capability_tests_20261001.json) pins that already-inspected log;
no sibling failure cause is inferred and no remote job is restarted.

Require positive native RSS only where `/proc/self/statm` exists. Keep the actual
pipeline export test on supported Unix platforms: when `/proc` is absent, it
must retain `rss_measurement="unavailable"`, zero RSS sentinels, stage timing,
counts and metadata. This does not accept zero as measured memory. Two new
tests confirm that the unchanged OrthoBench resource readers reject this
unavailable convention even when other numeric fields are valid.

Additional tests exercise disabled observation, page-to-byte conversion,
malformed and inaccessible reads, child token parsing, duplicate/cyclic process
traversal, failure evidence, propagated exceptions, SystemExit status and stopped
monitor threads. The explicit absent-proc fixture changes only `/proc/` reads
and existence checks; ordinary temporary-file I/O remains real. No global
autouse fixture, production fallback or weakened resource guard is introduced.

## Executed Validation

**54 focused Linux cases pass in 1.09 seconds**, zero failures/errors/skips,
including 17 new cases: metrics 17, retained-resource audit 10, factorial
assembler 27. The actual positive native RSS test passes on this host.

A separate process deliberately hides `/proc/` before test-module import.
The same 54-case scope gives **53 passed and one native RSS skip in 1.06
seconds**, no failures/errors. The pipeline export and both unavailable-memory
rejection cases still execute. This is a Linux capability simulation, not actual
macOS execution; the panels overlap and must not be added together.

Both JUnit reports, test-source pins and the three unchanged production-source
pins are retained in the receipt. Production files are byte-equal to the base
commit. The private environment is reused without installing or upgrading any
package. No scientific workload, service change, host contention poll, DGX
access or unrelated job action occurs.

At 02:44:06 UTC on 2 October, source-9674 CI is terminal failure in all five
test jobs, with CPU-wheel/docs success; no source-9674 log was inspected here.
Source-8b48 CI has CPU-wheel/docs success and five live test jobs. That run does
not contain this RSS correction. The prospective Bash-runtime correction and
this test change both still require actual remote confirmation.

## Remaining Boundaries

Sampled process-tree RSS is not an exact cgroup peak or unique physical memory.
Individual failed proc reads retain existing zero sentinels; the collector's
availability label depends on file existence, not proof that every descendant
was readable. These tests do not certify complete sampling or make historical
measurements controlled. Native benchmark admission, environmental handoff,
comparable timing, source/rights, remaining uncertainty and the versioned release
remain separate unfinished requirements. The publication goal stays active.
