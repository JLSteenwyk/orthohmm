# Public and Raw Benchmark Test Profiles

Public CI cannot provision the private SwissTrees source archives. It now uses
an explicitly named public profile rather than repeatedly presenting four
missing-input failures as the full test workflow. The original default test
and coverage targets still include those four cases. Their scientific,
source-validation, numerical and no-overwrite assertions are unchanged.

This is a test execution boundary, not completion of the publication release.
No benchmark scores, default parameters, raw inputs or native experiments change.

## Selection and Execution

The `raw_benchmark` marker labels exactly these four existing cases:

```text
tests/test_export_swiss_duplication_strata.py::test_export_and_no_overwrite
tests/test_export_swiss_duplication_strata.py::test_independent_rational_reproduction_of_export
tests/unit/test_swiss_strata_count_selection.py::test_explicit_counts_preserve_previous_rows_and_reject_wrong_hash[fragment]
tests/unit/test_swiss_strata_count_selection.py::test_explicit_counts_preserve_previous_rows_and_reject_wrong_hash[duplication]
```

No automatic missing-file skip, source replacement, mock validation or default
marker filtering is introduced. The descriptive/identity parameter cases,
configuration tests and synthetic exporter checks remain in the public profile.
An ordinary default collection still selects all four raw-data cases.

Before each macOS test job, the structured collection audit requires exactly
these four marked/deselected node IDs and every other collected node to remain
selected. Additional exclusions, changed markers, duplicates, collection errors
and collection-level skips fail the gate. The audit executes no test bodies
and does not establish that the selected tests pass. Its saved JSON includes
the complete selected inventory, collection diagnostics and interpreter identity.

The job names distinguish `public-fast` (also excludes `slow`) from
`public-coverage` (includes `slow`, excludes the four raw-data cases). Both
retain selection JSON and separate unit/integration JUnit receipts, including
failed outcomes. Existing execution-time platform/data-dependent skips remain
possible and are recorded in JUnit; a public-profile pass is not execution of
all biological/native regression cases or the complete raw-benchmark gate.
The Linux diagnostics and installed CPU-wheel jobs are unchanged.

With application/test dependencies and native test utilities available:

```bash
python -m benchmark_tools.audit_public_test_profile \
  --root "$PWD" --output /new/receipts/public-test-selection.json
make test.public.fast TEST_RECEIPT_DIR=/new/receipts
# Includes slow tests and emits coverage plus JUnit:
make test.public.coverage TEST_RECEIPT_DIR=/new/receipts
```

Run from the checkout root and use a fresh audit output path. The collection
inventory has no slow filtering; the fast execution command does. The audit
does not verify every runtime/transitive dependency, clear data rights or
replace the separate native-installation and scientific-admission checks.

## Full Raw Benchmark Gate

Keep using the [explicit checksum-bound raw-input handoff](SWISS_RAW_PYTEST_HANDOFF_20261002.md)
and [private archive restoration](SWISS_RAW_ARCHIVE_RESTORATION_20261002.md).
The unchanged default Make targets now accept `PYTEST_ARGS` to forward those
options, without adding a marker exclusion:

```bash
make test PYTEST_ARGS='--swiss-duplication-bindings /restored/duplication-inputs/bindings.json --swiss-duplication-bindings-sha256 TRUSTED_DUPLICATION_SHA256 --swiss-fragment-bindings /restored/fragment-inputs/bindings.json --swiss-fragment-bindings-sha256 TRUSTED_FRAGMENT_SHA256'
```

The uppercase digest tokens are placeholders for independently retained
digests, not hashes to calculate from untrusted received manifests. A real
run requires the complete legitimately supplied raw directories and compatible
test/native runtimes. No raw files or new redistribution permissions are
provided by this patch. Default `make test`, `make test.fast` and
`make test.coverage` do not omit the raw cases; only the explicitly public
targets do. Public targets do not consume `PYTEST_ARGS`.

## Executed Validation

The final combined five-module panel in the previously prepared private
Python 3.12.3 / pytest 9.1.1 environment passes **72 cases in 8.63 seconds**,
with zero errors, failures or skips. It includes 25 profile/audit cases, 13
binding-configuration cases, 20 dependency/Make cases and all 14 existing
SwissTrees exporter/count-selection cases. All four private cases execute
with the real restored raw inputs; no production validation is mocked.

The same scoped public selection passes **68 cases with four deselected**
in 3.02 seconds. These panels overlap and must not be added. Collection-only
verification of the complete checkout records **14,295 nodes: 14,291 selected,
four deselected, zero collection skips/errors**. That is not a full-suite pass.
Fresh subprocess checks also verify unfiltered/public/raw-only selection on
the two actual exporter modules without executing their tests.

Retain the first shared-interpreter failure: 53 passes/four failures in 4.05
seconds, caused by older Biopython/Pillow/setuptools and absent sqlglot. No
shared packages are upgraded. Its initial collection audit mistakenly allowed
two skipped sqlglot modules and omitted 13 cases. Preserve that output, add
explicit collection-skip rejection, and rerun the audit: the shared environment
now fails with both omitted modules identified, while the prepared private
environment succeeds. Twenty signal-parameter node labels also differ between
Python versions; these are renamings, not additional omitted cases.

The final combined JUnit is 11,108 bytes, SHA256
`e90e4898e552c044af96f4ccef04eb728231cb96d97dad55c56a918a37635728`.
The [machine-readable receipt](public_test_profile_review_20261002.json)
pins source, positive/negative collection evidence and all retained test attempts.

## Preceding CI and Remaining Work

At 12:55:49 UTC, exact-source `61082b22` run 37008346694 is terminal:
docs/Linux/wheel succeed and all five macOS test jobs fail. Its actual 3.13
job log is downloaded once at 12:57:17 UTC and verifies **14,142 passes/four
failures/119 skips/30 warnings in 386.33 seconds**. All 18 new accuracy-overview
cases pass; the four failures are the declared raw-data cases. Do not infer
sibling counts or successful full CI. This preceding run did not use the new
public profile; the amended workflow needs separate post-push confirmation.

Timing remains deferred without host polls/questions, DGX access or unrelated
job/service changes. Controlled resources, other-QfO uncertainty, transitive
runtime/data rights, complete executable release, archival deposition and final
manuscript/package reconciliation remain open. Neither a public-profile pass
nor the 72 focused cases establishes publication readiness.
