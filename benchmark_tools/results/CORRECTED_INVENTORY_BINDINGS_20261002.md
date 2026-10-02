# Corrected Inventory Test Bindings

The retained SwissTrees descriptor-inventory test followed four historical
absolute checkout paths. Bind only their in-memory test records to committed
files in the active checkout. Preserve every retained hash and byte count,
the complete descriptor/bin assertions, and production `check()` validation.
Neither scientific code, historical reports nor benchmark scores change.
Twenty negative cases cover four records with wrong hashes, wrong sizes,
same-size changes, truncation and absence. All must still fail validation.

## Verification

The first pre-change copied test passes by reading the original checkout;
it is not portability evidence. A second pre-change copied child, with
original-checkout Python opens forbidden, fails at the protocol path.
After the fix, all 44 affected/fixture cases pass in a copied isolated
Python 3.10.13 child in 0.61s: 1,816 staged files, 10,957,776 bytes and
15 staged project module origins. The original-checkout canary is blocked,
later original open events are zero, and child subprocesses are forbidden.
Temporary staging is removed. This Python-event guard is not OS containment,
native biological restoration or raw-input validation.

The complete local panel has **122 passes in 4.81s**, zero errors, failures
or skips: 37 inventory, seven record-fixture, 49 arithmetic-checker and
29 component-bundler cases. The earlier 44-case pass overlaps this panel.
An initial panel command named two nonexistent modules: exit 4, zero tests;
its receipt is retained, not counted as a successful check. Only 44 cases,
not all 122, were run in the new guarded copy.

```bash
python -m pytest -q tests/unit/test_corrected_swiss_sequence_strata.py tests/unit/test_retained_record_at_path.py tests/unit/test_check_swiss_descriptive_tables.py tests/unit/test_bundle_swiss_descriptive_tables.py
```

[Machine-readable identities and receipts](corrected_inventory_bindings_20261002.json)
pin the changed test, unchanged inputs/helper and successful/failed evidence.
The derived audit JSON is checked as a retained file; this does not
reconstruct its original FASTAs or independently verify raw biology.

## Remote And Raw Boundaries

Preceding source `6ea74c0d54727a7e41915b1154e75ccfe2a304f4`, run
36989749872, is terminal: all five macOS test jobs fail while Linux
diagnostics, wheel and docs succeed. Actual Python 3.13 job 110782933929
verifies that checkout and reports 13,943 passes, six failures, 118 skips,
zero errors and 30 warnings in 450.93s. All 21 frozen-lineage replay cases
pass, confirming that preceding fix. This inventory correction is not in
that source; do not infer new-patch remote confirmation or full CI success.

The remaining failures comprise this inventory metadata case, four raw
Swiss exporter/count-selection cases and one platform-specific live host
probe. Raw duplication admission lists 11 records/93,699,174 bytes: seven
main-Git tracked and four external/non-main-tracked records/93,650,891 bytes,
including the native Darwin image and QfO reference files. Fragment admission
lists 1,774 records/607,096,715 bytes: three main-Git tracked and 1,771
non-main-tracked records/606,389,322 bytes. These dependencies exist locally
but were only inventoried, not re-hashed or re-admitted here. Normal checkout
portability requires more than rebinding the first failing absolute path.
Do not replace raw-source gates with derived replay or skip them to green CI.

Timing remains deferred. The user's quiet-window answer imposes no immediate
coordination requirement: continue correctness/reproducibility work without
host contention polling, DGX access or unrelated process/service actions.
Controlled resource evidence, remaining uncertainty/provenance, rights,
complete release and archival deposition remain open. The publication goal
is active and incomplete.
