# Retained Overhead Replay Metadata

Four archival replay tests compare current file records with historical
absolute checkout paths. They pass in the original checkout, but all four
fail after copying the same source and data into a different directory.
The failures concern path metadata, not changed file bytes or scientific work.
This correction makes those tests meaningful outside the original checkout.

## Explicit Test Bindings

Use the existing, non-autouse `retained_record_at_path` fixture on only the
expected in-memory source/input records stamped by each unchanged builder:
three frontier records, two lineage records, two pressure records, and the
quiet CPU diagnostic's source, input and nine helper records. Do not rewrite
archived files, normalize arbitrary paths, monkeypatch production recording,
or repin hashes. Embedded historical commands, input paths, method order,
resource settings, counters, screening flags and admission decisions retain
complete nested equality. Pressure-plan disk digest checking is unchanged.

Strengthen the fixture's mutation tests: same-size changed bytes, truncation,
wrong byte length and wrong SHA-256 all remain rejected, alongside existing
changed/missing cases. The initial same-size test text was one byte short;
the failed 127-pass/one-failure receipt is retained. Change exactly one byte
before the successful verification. No production source changes.

## Verification

The four original-checkout tests pass before modification; that is not a
failed control. A fresh copied-tree pre-change run fails all four. After
correction, **128 local tests pass in 10.65s**, with zero failures, errors or
skips. This includes the five affected modules and the existing frontier
counter/replay regressions. Earlier runs overlap, not independent replications.

A new Python 3.10 isolated child passes **all 69 cases in the five affected
modules in 1.02s**. It uses 1,825 copied files/11,294,417 bytes, including
unchanged source and retained inputs. All 37 loaded project module origins
are staged. A Python audit guard rejects an original-checkout read canary,
observes zero subsequent original-checkout open events and forbids child
subprocesses. Temporary staging is removed. This is Python-event guarding,
not OS containment, cross-host execution, all-suite testing or full native
scientific restoration. Shared installed test dependencies remain external.

[Source, input, test and receipt identities](retained_overhead_metadata_20261002.json)
preserve the successful checks and both kinds of failure. Generated receipts
remain in ignored local work; the committed compact record pins their bytes.

## Preceding Remote Evidence

At 08:32:46 UTC, automatic source-0ab6f940 run 36983078805 is terminal:
all five macOS test jobs fail, while Linux diagnostics, CPU wheel and docs
succeed. Its actual Python 3.13 fast log has 13,904 passes, 14 failures,
118 skips, zero errors and 30 warnings in 392.84s. The four metadata cases
still fail; this new correction is not in that source. The 22 frontier
counter and 17 frontier-native cases pass, and the earlier portable exact
replay failures are absent. This confirms that bounded preceding correction,
not all 172 local cases remotely, complete CI or the new metadata patch.

Observe the next automatic CI run after committing; do not restart or
resubmit completed handles. Remaining raw/provenance/platform failures,
other-QfO uncertainty, rights, controlled resource evidence, complete release
and public deposition stay open. Scientific implementation/settings/scores,
historical archives and production admission are unchanged. Timing stays
deferred, with no host-contention poll, scheduling question, DGX access,
unrelated process/service action or new biological computation.
