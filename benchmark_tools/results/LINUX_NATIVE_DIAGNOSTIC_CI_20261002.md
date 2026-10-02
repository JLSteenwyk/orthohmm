# Linux Native Diagnostic Coverage

## Failure And Coverage

The retained source-3485 macOS Python 3.13 fast log has six fresh graph
worker failures at missing `os.sched_getaffinity`, before graph observations.
Two high-CPM child tests fail without exposing captured stderr; they invoke
the same Linux-only worker, but that is not direct traceback confirmation.
Those diagnostics also require `/proc/self/maps` and verified affinity setting.
These are Linux execution requirements, not observed constructor mismatches.
The existing repeat-worker tests already declare the same capability gate.

Apply that explicit gate to the eight constructor/direct/high-CPM worker
cases only. Leave their pure validation tests running on every platform.
Do not emulate affinity, omit loaded-library evidence or relax production
diagnostic admission. Add a separate Ubuntu/Python 3.12 CI job running all
five graph/parser test modules. Its JUnit gate requires all eight constructor
and twelve repeat-worker cases and rejects any selected skip, failure or error.
Always retain the receipt, including failed runs. This is execution on a
capable runner, not conversion of unavailable native evidence into success.

The parser orchestrator's four mocked tests also attempt to hash a Linux
libc file during fixture setup on macOS. Bind only that exact record request
to a synthetic, hashed temporary file within these intercepted-child tests.
Preserve the logical path and all other real record calls. Add three cases
showing same-size mutation, truncation and deletion prevent child launch and
output creation. No real libc is substituted or loaded. Production parser
validation still requires its exact retained interpreter/libc records.

## Validation

**317 ordered local cases pass in 33.53s**, zero failures/errors/skips:
the 76-case diagnostic scope followed by the preceding 241-case portability
scope. Twenty native worker cases execute, with unchanged endpoint/weight,
isolate, import-mode, affinity and no-unplanned-optimizer assertions. This
includes twelve repeat-worker variants, not twenty-four. These are tiny
synthetic graphs; they do not rerun QfO or demonstrate full-scale stability.
The initial 76-case pass overlaps and is not an independent replication.

The first local exercise of the new JUnit gate incorrectly expected 24
repeat cases and rejected the valid 76-case report. Inspection established
the actual 2-by-2-by-3 matrix of twelve; correct the prospective gate before
commit. The corrected code accepts the real report and rejects five generated
negative variants: missing constructor, missing repeat, skip, failure and
error. Temporary generated copies are removed. For each of three decorators,
evaluate the unchanged capability expression under four shapes: available,
missing maps, missing getter and missing setter. All twelve checks pass.
YAML parses using local PyYAML; no new runtime dependency is added.

[Machine-readable identities and commands](linux_native_diagnostic_ci_20261002.json)
pin the workflow, tests, reports and already-retained failure log. Production
diagnostics/core remain unchanged relative to source-c7b. Historical scientific
receipts and archives are not repinned. The synthetic parser binding is not
native parser replay, and local Linux success is not new-runner confirmation.

At 07:37:05 UTC source-c7b run 36978663641 still has five live test jobs,
wheel/docs successes. Do not restart or infer results from that incomplete
observation. This new workflow is not in that source. Actual execution of the
new Linux job and new macOS gate remains to be checked after push.

At 07:40:38 UTC the same preceding run is terminal: five test failures,
wheel/docs successes. Its actual Python 3.13 fast log confirms all 241 prior
archive/allocation/checker/cleanup/FAS/provider cases pass, including the two
previous fixture fixes. Overall: 13,890 pass, 28 fail, zero errors, 110 skip,
30 warnings, 481.23s. This confirms that prior scope on macOS, not the new
Linux job, full suite or sibling jobs. Preserve the log and do not restart.

## Remaining Work

Complete native workflow restoration, workstation/raw/platform failures,
dependency/source and data rights, other-QfO uncertainty, controlled resources
and versioned public deposition remain open. This work neither improves nor
changes scientific settings/scores. No biological inference, bootstrap,
annotation, scoring, raw search or old archive regeneration occurs. Timing
stays deferred, with no host-contention poll, scheduling question, DGX access
or unrelated process/service action. Publication readiness is not established.
