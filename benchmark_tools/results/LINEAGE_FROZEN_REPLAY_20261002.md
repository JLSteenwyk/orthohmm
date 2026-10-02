# Frozen Lineage Control Replay

Historical read-crossing controls must be replayed with the exact source
that produced them, not current helpers under old source pins. All eight
deployed source records match commit `6599c6e`; only the current frontier
helper differs. In the original checkout, the complete replay consequently
fails its local dependency guard. On macOS, the historical scheduler command
path fails first and masks the intended source-negative test.

## Test Staging

Keep production replay and every validation gate unchanged. The test fixture
creates a temporary shared Git clone without checkout, materializes eleven
historical files from the pinned commit, and copies the unchanged current
replay driver. The eight deployed source pins remain authoritative; three
additional files provide package, recording-helper and batch-path context.
The clone intentionally shares original Git objects and requires retained
history. It is not a self-contained archival restoration.

Run complete/source/scheduler checks in a fresh isolated interpreter importing
the staged dependencies. Bind only the exact single `Command` line in a
temporary scheduler copy to that checkout. Preserve all other scheduler bytes,
allocation/status fields, the 35 original raw files, archived service commands,
scope identities, source digests and complete control results. No production
record, historical receipt or scientific command is rewritten or repinned.

Add six negative cases: wrong scheduler command, same-size/truncated/missing
protocol bytes, same-size imported Python helper changes, and substitution of
the current frontier helper. Existing wrong source hash/allocation, incomplete
archive and raw-trial mutations remain failures. Expose child stderr on failure.

## Verification

Before change, one selected local case passes and complete replay fails its
source guard. The local scheduler path coincides, so this is not the same
first failing gate as macOS. The initial 19-case and 132-case successful panels
overlap the final check. **134 focused tests pass in 3.04s**, zero failures,
errors or skips, including all 21 lineage cases and nine fresh-child cases.
The panel also preserves the preceding scaling-command/source tests.

An actual fresh system Python 3.12.3 child exactly reproduces all three nested
control results, preserving signed residuals and false scientific/environmental
admission flags. Its 36 evidence records comprise 35 raw files plus the
temporarily rebound scheduler. All eleven loaded project module origins are
staged. A Python open-event guard rejects the original-checkout canary and
observes zero later original Python open events; only frozen-source `git show`
child commands are permitted. Temporary external staging is removed.

Git itself still reads shared original objects outside Python's open-event
observation. This is not OS containment, original-repository independence,
cross-host execution, full native scientific restoration or controlled timing.
No live service, CPU burn, scheduler submission or host-counter probe runs.
[Source, raw inventory, test, receipt and log pins](lineage_frozen_replay_20261002.json)
retain the failure and bounded successful evidence.

Run the test from a checkout with the pinned historical commit available:

```bash
python -m pytest -q tests/unit/test_replay_lineage_read_crossing_control.py
```

## Preceding Remote Evidence

At 09:17:30 UTC source-3bea410f automatic run 36986995236 is terminal: all
five macOS test jobs fail, while Linux diagnostics/wheel/docs succeed. Its
actual Python 3.13 log verifies that checkout SHA and reports 13,935 passes,
eight failures, 118 skips, zero errors and 30 warnings in 364.21s. All 120
preceding command/metadata cases pass. The two lineage failures remain;
this new staging patch is not in that source. Do not infer new-patch remote
confirmation, other interpreter outcomes or full CI from the bounded log.

Scientific implementation/settings/scores, production code and archived raw
evidence stay unchanged relative to 3bea410f. Timing remains deferred without
contention questions/polls, DGX access or unrelated process/service actions.
Remaining raw/platform/provenance failures, other-QfO uncertainty, rights,
controlled resource evidence, complete release and public deposition remain
open. The original publication goal stays active and incomplete.
