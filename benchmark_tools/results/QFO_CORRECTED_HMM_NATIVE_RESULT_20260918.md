# Corrected QfO High-Sensitivity Native Result

## Admitted Scope

Primary inference `21706_0` (raw job `21707`) completed with exit `0:0`
after scheduler elapsed time `19:41:39`, using 32 allocated CPUs on `bizon`.
Native admission `21720` and replay preparation `21722` also completed
with exit `0:0`. This admits native evidence for downstream analysis,
**not biological accuracy, search completeness, or controlled timing**.

The retained [admission report](qfo_corrected_high_sensitivity_admission_21720.json)
has SHA-256
`845874021057da54cd61c20763b00280c1a5569c5f451676e52f16043dacccff`.
It binds the frozen command, runtime availability checks, source, corrected
inputs, native output inventory, metrics and numeric checkpoint. A fresh
rehash verified all 107 distinct records referenced by admission/content
checks and replay preparation.
Rerunning the frozen admission script against raw job `21707` reproduced
the entire report byte-for-byte (`cmp` exit `0`). The independent rerun is
retained at
`benchmarks/work/qfo_corrected_high_sensitivity_admission_recheck_20260918.json`.

| Validated native property | Value |
| --- | ---: |
| Species | 78 |
| Proteins covered exactly once | 984,137 |
| Final groups | 391,908 |
| Groups with one member | 302,763 |
| Singletons added to complete the raw partition | 0 |
| Numeric checkpoint hits | 90,687,327 |
| Checkpoint self-hits | 983,835 |

One-member groups are counted directly from the partition. They must not
be confused with the zero *added* singleton count. Checkpoint integrity
does not prove hit completeness or independently validate HMM scores.

Native group-file SHA-256:
`eaace6c6d0feb6bf441c1bcac8c9029c208dfbc34541fcad781ce13afa44196b`.
Checkpoint-manifest SHA-256:
`30d1718e0ad1ec080b9783b0f94f3e8855063459c4b2999d86f93eb8d1b4c8e1`.

## Resource Evidence

Native metrics record wall time `70886.239038` seconds and sampled summed
process-tree peak RSS `18245820416` bytes (`16.992744` GiB). This is Linux
process-tree RSS sampling, not an exclusive cgroup memory measurement.
It can miss unsampled peaks and double-count shared pages. The metrics
file was hashed after completion; it is a sibling of the inventoried native
output directory. Its SHA-256 is
`ba0f65d24ad594c802e3e38ddad04637f8f88694a3e8371042c17274d2d343e5`.
The run used the shared host; these observations do not establish a speed
or memory advantage over another method. Scheduler elapsed time includes
wrapper work and is distinct from native metrics wall time.

## Downstream State

The retained [replay preparation](qfo_corrected_replay_preparation_21722.json)
has SHA-256
`e7657732bc94438fb0602f3246680308a075b9fda1c68909ddfac1138b164c0c`.
Its `corrected_replay_command_frozen_unrun` status and
`execution_authorized: false` describe preparation: the native command
must not run without the separately frozen checked dispatcher. They are
not live scheduler status. Checked dispatcher `21756` is now running,
as is hit-coverage analysis `21793`. OrthoFinder `21706_1` has started.

Replay must check all four clustering boundaries and compare its final
partition with this native result before candidate expansion and the eight
P/C/R cells are admitted. Conversion, scoring and paired uncertainty remain
separate gates. No historical score is transferred to the corrected input,
and no new HMM accuracy score is claimed here. The scientific core remains
the frozen `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` revision.
