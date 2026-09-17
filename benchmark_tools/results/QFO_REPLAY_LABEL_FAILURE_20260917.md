# QfO Replay Wrapper Label Failure

Job21333 ended FAILED1:0 after51:17. This does not mean independent scientific
admission succeeded. The wrapper raised `Incomplete full replay or missing
profile construction` after its child replay returned exit0.

Inspection identifies a concrete label-contract mismatch:
`run_qfo_checked_full_replay.py` expects `profiles` and `profiles_refined`, but
the frozen scientific `replay_high_sensitivity.py:326` emits `strict_profiles`
and `strict_profiles_refined`. Its other two labels are `multipass` and
`multipass_refined`. The independent admission script currently repeats the same
incorrect label expectation and requires parent success, so it must not be
used unchanged to claim this failed job was completed successfully.

The preserved worker reports four checked, exit0 calls, each covering976504
genes: initial349898groups, multipass308190, profile_base334621,
profile_expanded305656. The returned replay reports57883profiles built,
19519955profile candidates,1966961significant profile hits and98085strict
profile edges. These are raw execution reports, not independently admitted
accuracy findings or proof of final coverage.

Evidence retained without alteration:

- [Failed parent](qfo_checked_full_replay_v2_label_failure_20260917.json), SHA256
  6a6b682e2e86566e59ddf8748f6eeb475eab976ee35c988df73bf1b75d4c100a.
- [Returned replay](qfo_checked_full_replay_v2_returned_20260917.json), SHA256
  ac0e69463a63fd6a4431ad6118bdc6a5835187d2b34c4f691ddfbdd3519d9885.
- Local `benchmarks/results/qfo_checked_full_replay_v2/checked_worker.json`,
  SHA256bed6806a82cadd6abcd828c0adca232040052580b31db9b0b997ba57010dca73.

Next: add a separately recorded, strictly scoped postflight recovery that
accepts only this exact failed-parent evidence, preserves the scheduler failure,
and validates the correct frozen label contract. Recheck all native payloads,
graph fingerprints, source/runtime/input identity, final coverage, initial repeat
agreement and historical partition disagreement. Do not overwrite the failed
report, weaken unrelated checks, choose partitions by accuracy, or rerun the
51-minute inference merely to repair bookkeeping. The original parent failed
before its after-input audit; a fresh audit must be labeled retrospective.

## Recovery Implementation

The corrected auditor recognizes the frozen `strict_profiles` labels, including
their partition-copy and historical-reference checks. An explicit
`--recover-label-failure` path accepts only versionv2 and the three exact SHA256
records above. It requires the original schedulerFAILED1:0 record, native
workerexit0, allfour checked clustering calls, the precise wrapper error and
absence of the postflight fields that the failed parent never produced.
Other failures are not recoverable through this path.

The auditor independently reconstructs the saved graphs, verifies all native
observations and sources, checks final partition coverage and copy identity,
performs a fresh retrospective input audit, and retains initial/historical
partition comparisons. It does not fabricate parent completion or modify old
outputs. On success its distinct status is
`checked_full_replay_recovered_verified`, with original failure and retrospective
limitations explicit. Normal admission still requires a successful parent and
completed scheduler record. The unexecuted current runner's label guard is also
corrected; frozen inference executors are unchanged.

Thirty-seven focused auditor tests pass, including exact recovery failure gates,
rejection of old labels, failed child calls and nonmatching scheduler states.
The pinned real reports pass the corrected inventory check only; full native
admission has not yet completed. Prepared four-CPU64GiB audit-only batch to
validate existing outputs without running clustering or profile construction.
