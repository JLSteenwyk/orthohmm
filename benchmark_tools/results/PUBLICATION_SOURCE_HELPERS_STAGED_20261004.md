# Source Helper Export Fix Staged

The [unapplied patch](publication_source_helpers_staged_20261004.patch) implements
the source-export repair in an isolated copy, preserving the live timing
recipe. It adds explicit `--include-result-helpers` selection and schema3;
historical default/profile exports and schema1/2 verification remain unchanged.
It exports committed top-level results-module Python sources, not JSON result
data, FASTA, nested helper directories or uncommitted working-tree bytes.
Scientific source and the separate setup-only overlay retain their revisions.

All187 focused legacy/new tests pass in35.00s. Fixture exports relocate without
the original checkout or Git, copied isolated CLI verification succeeds, helper
imports resolve through namespace/relative imports and an exported test collects.
The new tests cover all20 observed omissions and refuse ambiguous selection
flags or historical-schema reinterpretation. These are fixture tests, not
actual full-study source export, installed native execution or runtime/data
closure. Current code still uses the historical selector.

The initial run retains35 failures/160 passes: the new fixture accidentally
committed the old fixture's deliberately invalid working-tree probe. Correct
only the fixture's staging list and use the proper support-profile fixture;
the next187-case execution passes. The final tests are made self-contained,
without depending on the gap JSON, and the full187-case suite passes again.
All three JUnit files remain local/ignored and are pinned in the
[validation receipt](publication_source_helpers_staged_validation_20261004.json).

The patch passes `git apply --check` without changing the live project.
Application in another isolated directory reproduces both tested source files
byte-for-byte. Base source remains15,388 bytes,SHA-256
`b199863d32a2370d8ce0b20dfbeb2e082d3486a0fc89a0dafd1b568b2bb0cc32`.
Patch is14,069 bytes,SHA-256
`ce40376209667428bf6f7d947a71fb3cec12fd556538cc4221d8f04e67a71fb0`.

Integrate only after all27 frozen timing identities are terminal and reviewed,
rechecking the base before applying either file. Retain the tested new unit
file as `tests/unit/test_bundle_publication_source_helpers.py`; then commit
the actual helper before building a final committed-source export with the
explicit option. Validate relocated final workflow imports/collection and
reconcile final manuscript/reporting/archive evidence. This prepared fix
neither completes publication packaging nor authorizes timing retries.
