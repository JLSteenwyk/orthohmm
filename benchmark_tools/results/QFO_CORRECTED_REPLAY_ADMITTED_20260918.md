# Corrected QfO Replay Admission

Replay job 21756 completed with exit code zero in 50:28 on the shared local
host. Independent admission job 21757 completed with exit code zero in 1:17.
The [admission receipt](qfo_corrected_replay_admission_21757.json) checks 419
provenance records, the frozen runtime, scheduler identity, four native
clustering boundaries and complete partition coverage.

The final replay and corrected native high-sensitivity result have identical
partitions: 391,908 groups, zero native-only groups and zero replay-only
groups across the 984,137-gene universe. This is corrected-input evidence;
it does not retroactively establish equivalence for older historical runs.

| Retained output | Groups |
| --- | ---: |
| Multipass | 308,697 |
| Multipass refined | 394,328 |
| Strict profiles | 306,051 |
| Strict profiles refined | 391,908 |

These outputs are distinct from the four checked native clustering calls:
initial, multipass, profile base and profile expanded. Every native call has
an exit-zero checked execution receipt. Refinement outputs are checked
separately by the admission workflow.

## Evidence

The retained admission JSON is an unmodified copy of
`benchmarks/work/qfo_corrected_replay_admission_20260918.json` (108,715 bytes).
Both copies have SHA-256
`1eec5ffb675fff234e5ce0db7e65abfa9bec09bdea8f5bd13ca72ad72d7ec40e`.
The completed parent report at
`benchmarks/results/qfo_corrected_checked_replay_v1/results.json` has SHA-256
`c2217dacdaf1a6e184a21eed4a432535d08c4b81d77c2076d6a8901212857190`.
The final replay partition SHA-256 is
`d21bf01bae3bf2392818cc23babbc6fc1af0d6a0d1f144831ce05dfa9ec85e77`.

Admission independently rereads partitions and compares memberships rather
than trusting the parent's equality flag. It also rechecks bound file
hashes. No source executor, default parameter or output was changed to
obtain agreement, and no duplicate replay was launched.

## Remaining Work

Candidate preparation job 21758 started after successful admission; its
independent admission and downstream reconciliation/scoring remain pending.
The replay receipt explicitly sets `accuracy_evaluated=false` and
`publication_ready=false`. This milestone validates reusable intermediate
outputs, not biological accuracy, independent generalization or comparative
runtime. Cached shared-host execution is not dedicated end-to-end timing.
