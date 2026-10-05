# Native Factorial 22427: Scientific Recovery, Failed Wrapper

## Actual Outcome

Job **22427** is terminal **FAILED 1:0**, allocation elapsed **40:59**.
Its native wrapper also exited 1, without timeout. The scientific pipeline
completed all expected P0C0R0 stages and produced complete metrics and
partitions over **251,378 proteins from 12 OrthoBench species**.

The wrapper failed afterward while publishing its final execution receipt.
It called a create-once hard-link writer on the existing running receipt,
raising `FileExistsError`. The original running receipt and the complete
`.pending` receipt remain untouched, alongside the traceback and all raw
measurement files. This is not relabeled as a successful scheduler/native
command or a biological inference failure.

No search, clustering, tree construction or benchmark inference was repeated.
The broader publication goal remains active and incomplete.

## Independent Review And Recovery

The [terminal review](native_factorial_terminal_review_22427.json) binds the
original request/plan to fresh terminal accounting, independently checks the
runtime before/after records and actual terminal inventories, replays all raw
resource observations and re-evaluates typed process/pressure streams. It
also reproduces the launch-capacity, process-attribution and scheduler-budget
checks. Sampled shared-host evidence passes; no isolation claim is made.

The first terminal review, under Python 3.12, failed exact host-process replay.
At interval 1, the retained foreign-CPU sum is `41.37489194558179`, versus
`41.374891945581766` when recomputed in Python 3.12: difference
`-2.1316282072803006e-14`. Its failure and runtime review remain at
`benchmarks/work/native_factorial_terminal_review_22427`. No raw value was
rounded, changed or accepted with a newly chosen tolerance.

Use a separate [Python 3.10 review environment](native_factorial_review_interpreter_20261004.json),
matching the collector's interpreter generation, with NumPy 2.2.6,
Biopython 1.87 and psutil 7.2.2 installed from cached wheels. This changes
neither the frozen native nor private controller runtime. Initial missing-Bio
and missing-psutil import failures did not run inference or create terminal
review outputs. The matching interpreter reproduces the full raw review.
Retained work is at `benchmarks/work/native_factorial_terminal_review_22427_matched`.

A separately scoped [scientific recovery](native_factorial_recovery_22427.json)
checks the exact diagnosed executor/writer versions, initial and pending
receipt identities, completed metric/stage/count agreement and the final
traceback. It then independently [validates native outputs](native_factorial_outputs_22427.json)
and [rescores all 70 frozen RefOGs](native_factorial_score_22427.json), including
low-certainty exclusions. The original cached score is independently
reproduced before comparison. Recovery preserves the original failure.

## Scientific Result

This is the **P0C0R0 ablation**, not the main high-sensitivity or phylogenetic
method: profile expansion, candidate-family expansion and reconciliation are
off. Sensitive initial HMM search, k=4, candidate cap=100, multipass graph and
Leiden seed=4 remain unchanged.

| Endpoint | Full Native P0C0R0 | Retained Cached P0C0R0 |
|---|---:|---:|
| OrthoBench weighted F1 | 0.697633878435 | 0.697633878435 |
| Precision | 0.788685918309 | 0.788685918309 |
| Recall | 0.625429437185 | 0.625429437185 |
| Orthogroups | 63,245 | 63,245 |

The **entire gene partition is identical**, ignoring group labels: zero
changed groups/genes, not merely equal scores or equal counts. This provides
a full-native check of one earlier cached factorial cell. It does not prove
equivalence of other cells, paired superiority, independent generalization,
or publication readiness. OrthoBench is development-exposed evidence.

## Resources And Contention

The [replayed resource endpoints](native_factorial_resources_22427.json) are:

| Endpoint | Observed Value | Scope |
|---|---:|---|
| Native wall | 2,289.764741234 s | Complete native-wrapper command interval |
| CPU | 67,271.677476 s | Native-task subtree counter bracket, including wrapper |
| Peak memory | 4,954,263,552 bytes | Native-step lifetime peak, including launcher |

These are **failed-wrapper attempt** observations after complete scientific
inference, not clean-success timings. Allocation elapsed includes preparation
and reporting and is not substituted for native elapsed time. Internal
sampled process-tree RSS is a different memory scope and is not substituted
for the native-step peak. Whole-job CPU/peak remain unclaimed.

Whole-run monitoring observed up to **41.88460354100697** competing CPU-core
equivalents. Timings were collected on a shared Threadripper while other
analyses were running. CPU, memory-bandwidth and I/O contention may affect
elapsed times by an unknown, potentially tool-dependent amount. Matching
resource limits does not establish isolated tool performance. No definitive
speed ranking, causal component overhead or background correction is claimed.

## Narrow Receipt Fix And Next Work

After the job and independent reviews finished, replace only the wrapper's
terminal-receipt write with an atomic update of its owned running receipt.
Check the prior content, reject symlinks/foreign changes and unexpected
temporary-file collisions, flush the new JSON, then atomically replace the
receipt. The generic create-once writer remains unchanged for immutable
evidence. Scientific configuration, native core, adapter and resources are
unchanged. Preserve the [executed failing source](native_factorial_failed_source_22427/run_native_factorial_cost.py)
byte-for-byte; its SHA remains `07a926e38ba49144b56bf9a7d23478d3dd050fc3a14b98e1ce8889d1f98a8d20`.
Historical source references describe executed bytes, not current-source
identity after this repair; use that explicit archived copy for this one
changed path when validating historical evidence.

The [final test report](native_factorial_terminal_recovery_tests_20261004.xml)
records **496 passed**, zero failures/errors/skips, in 6.75 s. Tests cover
terminal joins, native failures/timeouts, runtime and launch/replay tampering,
all eight retained native diagnostic outputs, exact receipt-failure recovery
and successful/failed atomic receipt transitions. The preceding run retained
one historical-source test failure (495 passed); it was corrected to verify
the archived executed source, not falsely require the repaired source to
have the historical hash.

Next freeze a prospective execution/source amendment for the **next different
identity**, preserving this failed attempt and its recovered scientific result.
The earlier plan pins the failing executor and must not be reused as a valid
launch plan after the repair. Its completed index-0 review also belongs to
that historical plan: explicitly validate any cross-plan sequential-history
adoption rather than editing its plan/hash or silently restarting index 0.
Fresh capacity and accounting checks remain required, but neither contention
alone nor lack of a quiet window/DGX blocks progress. No new identity is
launched at this milestone; no unrelated analysis or service is changed.
