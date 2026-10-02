# Affinity Test Capability Boundaries

Make synthetic affinity guards independent of Linux API availability while
retaining genuine native execution checks on supported hosts. This is test
portability and clearer capability reporting, not a production fallback or
controlled-resource result. Continue from pushed `bec79e1d`.

## Observed CI

Inspect source-57a8391f macOS Python 3.11 job 110671849352 once. Its retained
log reports **13,489 passed, 159 failed, 28 errors, 96 skipped, 30 warnings in
462.83 seconds**. Measurement (21), saved-evidence replay (10), cgroup-frontier
(20), shared boot fixture (two), rendering (61) and count reproduction (seven)
cases pass. Fake-service capture still fails: `/bin/true` is unverified in the
synthetic reply. This does not diagnose why that host file was unverified.
Eight other modules have 49 affinity-related failures/errors, including native
child `preexec_fn` failures. Broad CI is not green, and sibling causes are not
inferred. All five source-57 test jobs are terminal failure at 02:15:54 UTC;
CPU-wheel/docs succeeded. Source-bec run 36954409877 has CPU-wheel/docs success
and five live test jobs at the subsequent recorded snapshot. Do not restart them.

## Scoped Corrections

Synthetic tests explicitly install their fake affinity APIs with `raising=False`;
fake panel launch records receive a fixed 20-CPU fixture, not the real machine's
mask. Negative resource-bound tests use a one-CPU fake mask. No autouse/global
affinity fallback is introduced. Guard expectations, source pins, allocation
limits, outcome retention and failed-run handling stay intact.

Genuine Linux workload, resource-limit child and loaded-library worker checks
require their actual capabilities. They skip explicitly when unavailable;
their native implementations are not mocked. Two new cases prove that an absent
affinity API still fails instead of synthesizing an allocation, and worker
creation is unreachable in that state. Fake-service replies now declare a
temporary regular file that is hashed but never launched, instead of `/bin/true`.
Real file-identity policy and secret-exclusion checks remain unchanged.

## Executed Validation

On this Linux host **203 cases pass in 16.08 seconds**, no failures/errors/skips,
including the 21 bounded native cases. Separately remove both affinity APIs
before importing test modules: **182 pass, 21 explicitly skip in 4.15 seconds**,
no failures/errors. The skip inventory is exactly 12 loaded-library worker,
seven native workload and two resource-limit child cases. Synthetic guards,
missing-API rejection and fake-service identity still execute. This is an
API-absence simulation on Linux, not execution on macOS. The panels overlap
and are not additive. All nine production modules remain byte-equal to the
prior commit. [Receipt and source/JUnit pins](ci_affinity_fixtures_20261001.json).

Commit/push, then observe the actual new automatic CI without manual retries.
New remote corrections remain unconfirmed. Historical input/path failures,
broader CI and final executable release portability remain separate work.
No DGX connection, host-contention poll/question, unrelated job/service action,
shared-package upgrade, native benchmark rerun or scientific/default/score change
occurs here. Tiny owned fixtures are not calibration or controlled timing.
The full publication goal remains active: comparable resources, source/rights,
uncertainty gaps, release and archival deposition remain unfinished.
