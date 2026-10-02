# Python 3.13 Test Fixture Compatibility

## Observed Failures

The retained source-055 macOS Python 3.13 fast log shows two test-fixture
problems, not demonstrated production failures. The oversized archive fixture
calls `TarFile.addfile` without a body for a nonempty regular member. Python
3.13 rejects that producer call before OrthoHMM's consumer executes. The
[official tarfile documentation](https://docs.python.org/3.13/library/tarfile.html#tarfile.TarFile.addfile)
explicitly records this version change.

Eight allocation traces use `set <= frame.f_locals.keys()`; the actual log
shows a list-form key collection and TypeError. [PEP 667](https://peps.python.org/pep-0667/)
documents the new frame-locals proxy; the specific key-shape evidence is the
CI traceback, not an assumption that every mapping has a set-like keys view.

## Test-Only Correction

Construct valid tar members with their complete declared bodies, including the
4-MiB-plus-one-byte negative case. Retain the production 4-MiB budget unchanged.
Require the specific consumer rejection for each of six existing faults and
forbid arithmetic execution. Add two real-body cases: a single member exactly
at the budget reaches missing-inventory rejection; two individually valid
members exceed the cumulative budget by two bytes and fail the budget guard.
No output is accepted. Compressible synthetic bodies keep fixtures bounded;
no raw data or large artifact is committed.

Use `fields.issubset(local_keys)` for the memory trace. Run each of the same
eight native allocation seeds with both runtime-native keys and explicitly
list-form keys. Retain the exact named-array `.nbytes`, snapshot bounds, input
logical-byte assertions and trace restoration. This is not a total peak-RAM
bound or graph-feasibility admission; the production estimator and accuracy
implementation are unchanged.

## Validation

**241 local tests pass in 13.35s**, zero failures/errors/skips: 28 memory,
29 bundler, 49 checker, 15 cleanup-helper, 37 measurement, eight WGD, ten
measurement-audit, 33 FAS and 32 provider cases. The intermediate 57-case
pass overlaps. Before replacing the set operator, the expanded memory panel
has 20 passes and eight failures, all eight forced-list cases reproducing
the same TypeError. Preserve that failed report. The sixteen allocation
cases are eight seeds tested twice, not sixteen independent datasets.

Local execution uses Python 3.10.13. No local Python 3.13 interpreter was
available or installed. Valid tar construction and list-key simulation are
focused compatibility evidence, not actual execution of this new patch on
Python 3.13. [Machine-readable pins and commands](py313_fixture_compatibility_20261002.json)
bind test sources, reports and the actual remote logs.

At 07:22:38 UTC the preceding source-3485 automatic run 36977059602 is
terminal: five test jobs fail, wheel/docs succeed. Inspect only its actual
macOS Python 3.13 fast log: all 135 previous cleanup/FAS/provider cases pass,
including the native TERM-ignoring descendant check. Overall: 13,871 passes,
37 failures, zero errors, 110 skips and 30 warnings in 485.96s. The same
one tar-producer and eight allocation-key failures remain; this new patch
is not in that source. No sibling/full success or exact kernel diagnosis
is inferred, and no CI job is restarted.

## Boundaries

Only two test files change. Scientific implementation/settings/scores,
production bundler/checker/estimator/helper and all historical archives and
receipts remain unchanged relative to source-3485. New-patch remote
confirmation and other workstation/raw/platform CI failures remain open.
No inference, bootstrap, scoring, annotation, archive regeneration or source
search occurs. Timing remains deferred, without a contention poll, scheduling
question, DGX access or unrelated workload/service action. Other-QfO
uncertainty, rights, comparable resources, full executable release and public
deposition remain unresolved; publication readiness is not established.
