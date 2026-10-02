# Remote Regression Confirmation

[CI run 37064634051](https://github.com/JLSteenwyk/orthohmm/actions/runs/37064634051)
is terminal and successful at exact source
`9c81953db35edb780f8546cf81a9133d8df4d129`. The final API observation records
eight completed successful jobs: four Python-version fast jobs, coverage,
Linux native diagnostics, CPU-wheel verification and docs. No jobs were
restarted. This is the glibc/ABI/assembly milestone, before citation metadata.

## Downloaded Execution Evidence

Download only the public-coverage artifact, ID 11252359194, 516,576 bytes.
Its SHA256 matches GitHub's published digest:
`a6200452af5de8cfa64ee635299fef6ccc2d8a59c0511a2e8a916de765b77672`.
Read JSON/XML directly from the ZIP; no archive extraction or test rerun.

| Execution receipt | Passed | Skipped | Failed/errors | Suite time |
| --- | ---: | ---: | ---: | ---: |
| Python 3.11.9 public unit | 14,878 | 119 | 0/0 | 592.330 s |
| Public integration | 5 | 0 | 0/0 | 202.898 s |
| Combined | 14,883 | 119 | 0/0 | Not a benchmark duration |

Selection records 15,006 collected nodes and exactly four declared private
cases deselected, with zero collection errors/skips/problems. All 15,002
selected IDs reconcile exactly to unique JUnit cases across both receipts.
Skipped executions remain accounted for, not reported as passes.

Five recent modules execute without skips: glibc guard 35, integrated workflow
31, ABI inventory 4, wheel ELF inventory 32 and runtime assembly 29 (131 total).
All ten corresponding source/test files currently match that Git snapshot
byte-for-byte. This does not compare every transitive dependency or claim that
macOS mock/fixture cases are actual private-controller Linux preflight.

First readback fails before output because splitting pytest IDs on every `::`
misinterprets 13 parameter labels. Correct the mapping using known module
paths and structured JUnit classname/name fields, preserving parameter text
literally. The receipts and source tests remain unchanged. Final readback
passes; its retained driver and the failure are identified in the
[machine-readable receipt](ci_glibc_confirmation_20261002.json).
Receipt: 8,294 bytes; SHA256
`e59df14415f960982efac40296e92d574175743ecee0cbd83ec7fbf387feef44`.
Execution counts are verified only for this downloaded artifact; sibling job
success does not establish sibling counts. No temporary download URLs are persisted.

## Remaining Publication Work

The public profile excludes four raw-source tests and retains 119 execution
skips; this is not the full raw/native regression gate or scientific admission.
The later citation commit and this documentation are outside that run's scope.
The original publication goal remains incomplete.

The [prospectively adopted resource scopes](THREADRIPPER_RESOURCE_SCOPE_DECISION_20260929.md)
already designate native inference wall time, bracket CPU and step-lifetime
peak memory as primary. Unavailable final whole-job teardown counters are
supplementary; they alone do not block those endpoints. Do not add scheduler
accounting changes as a prerequisite or label partial counters as final peaks.
Actual observer-overhead evaluation, real environmental handoff, current
source/runtime readiness and verified isolation still precede controlled runs.
No quiet-window question, host contention poll, timing/inference launch, DGX
access, service change or dated-archive rebuild occurs here.
