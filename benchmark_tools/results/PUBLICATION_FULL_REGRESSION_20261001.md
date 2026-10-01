# Complete Local Regression

At `c022ee78b100593f1e789b6842b3953186cbc5f7`, recorded in the [actual receipt](publication_full_refresh_20261001.json),
one full default-discovery invocation passed **13,660 tests**, with zero
failures, errors or skips. Pytest reported 22 warnings and 614.43 seconds;
the controller elapsed time was 622.271 seconds. These are shared-host test
durations, not controlled benchmark measurements.

The run includes all five integration cases, 232 top-level cases, the 13
SQL-parser cases previously hidden by two collection skips, and all ten
opt-in installed OrthoMCL fixture cases. Integration used explicit phmmer
and compatible MCL paths; the legacy OrthoMCL MCL was not substituted or
changed. The retained SQL parser reports version 30.20.0. No package was
installed or upgraded.

The receipt records the exact interpreter, command, scoped environment
overrides, package versions, start/finish times and process exit. Its local
JUnit and log remain at the recorded paths:

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| JUnit | 1,860,867 | `78c36d828796f0726de9d618a5bd1e7eb13d71e370a0ce7e4bd0f93d04fadc70` |
| Log | 15,692 | `5e0b970f8201e08ff03035967b313bd6276350b7bbe62ef20da6097000214ebf` |

The controller rechecked all 2,891 inventoried tracked source/test/fixture
files and selected executable bytes after the child terminated. Both
inventories were unchanged, the scoped tracked source status stayed clean,
and all 1,140 already-modified sample files retained their original observed
bytes. These inventories are pinned in the receipt. The committed receipt
is an exact copy of the actual controller result, not a reconstructed success
summary.

The warnings concern invalid escape sequences in frozen `parser.py` and
`writer.py`; they were not suppressed or repaired by altering the scientific
baseline. Native tests use small fixtures, not full benchmark datasets.
Passing this local run does not establish remote CI, complete dependency
closure, scientific accuracy, data rights, public deposition or publication
readiness. Completed analyses were not restarted, controlled timing remains
deferred, and no DGX, unrelated job/service or scientific default was changed.
