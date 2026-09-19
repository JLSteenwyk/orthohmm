# CPM QfO Assessment and Independent Admission

## Status

Submitted two serial assessment tasks (21982) and two serial independent
admission tasks (21984). All four are scheduler-confirmed PENDING on
dependencies. No new CPM score or uncertainty result is available or
admitted. The prespecified arms remain `cpm_low` (0.08) and `cpm_high`
(0.12), each with its own independently inferred phylogeny. This is
development-exposed parameter robustness, not a default-selection exercise.

## Frozen Assessment

- Executor commit: `43becaf4b2505a378edc5ad3111afcbf406e529b` (pushed).
- Worktree: `benchmarks/work/publication_qfo_cpm_assessment_v1`.
- Runner: `benchmark_tools/run_qfo_cpm_assessment.py`, SHA-256
  `8c5113e6e05d22c2259ce8aa5f22285aa077c4b5fdd7c85c2a5a71f26677a073`.
- Batch: `qfo_cpm_assessment_batch_20260919.sh`, SHA-256
  `0408e705d70046646d7d542fa428886b7a08881dab52a56b1918f74a4a68e008`.
- Allocation: 8 CPUs, 96 GiB, 24 hours, bizon, no requeue, array `0-1%1`.
- Dependency: `afterany:21978,aftercorr:21978`.

The runner first requires terminal successful pair conversion, with exact
raw job ID, task index, native-pair semantics and zero mapping loss. It
pins converter source/worktree, revalidates the CPM arm context and FASTAs,
checks the frozen QfO mapping/environment and all recorded inputs, and
uses the existing frozen six-endpoint command builder. Namespaces are
exclusive: `benchmarks/results/qfo_cpm_assessment_v1/<arm>`,
`qfo_benchmark/w/qcpv<index>` and
`qfo_benchmark/scoring/cpm_v1_<index>`. Existing paths, including broken
symlinks, are rejected. Preflight, logs and failure/success reports are
retained; successful process exit is not accuracy admission.

```bash
sbatch --parsable benchmark_tools/results/qfo_cpm_assessment_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_assessment_v1 \
  43becaf4b2505a378edc5ad3111afcbf406e529b
```

## Independent Admission

- Executor commit: `e4b1305b280b9fc6694b1071327cbc6e3c0442e0` (pushed).
- Worktree: `benchmarks/work/publication_qfo_cpm_assessment_admission_v1`.
- Validator: `benchmark_tools/admit_qfo_cpm_assessment.py`, SHA-256
  `907e5abbf7ec5b60781dfd48d95452bb67578662b1aefa44af8be131b7542c66`.
- Batch: `qfo_cpm_assessment_admission_batch_20260919.sh`, SHA-256
  `6d480f6fd3486ef42d2f03964f45407ac09423d44a3017f6206e8cf1e693be75`.
- Allocation: 2 CPUs, 64 GiB, four hours, bizon, no requeue, array `0-1%1`.
- Dependency: `afterany:21982,aftercorr:21982`.
- Reports: `benchmarks/work/qfo_cpm_assessment_admission_21984_<index>.json`.

The validator independently gates on terminal assessment success, verifies
preflight/postflight identity, converter and assessment checkout/source,
arm context, inputs, command, environment, exact output inventory and
unique native trace. Shared tested QfO validators check native task
completion, participant identity and all six endpoints' raw arithmetic.
Input/output hashes are rechecked before the admission report is written.
No benchmark label or parameter is changed. Full publication readiness and
paired uncertainty remain false or unmet even after eventual score admission.

```bash
sbatch --parsable benchmark_tools/results/qfo_cpm_assessment_admission_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_assessment_admission_v1 \
  e4b1305b280b9fc6694b1071327cbc6e3c0442e0
```

Whole-array terminal barriers preserve stable accounting snapshots while
corresponding-task success gates prevent a failed arm from entering scoring
or admission. Failure artifacts remain retained; there is no implicit retry.

## Verification and Remaining Work

Assessment/conversion regression suite: 115 tests passed in 0.67 seconds.
Admission/assessment/shared-endpoint suite: 142 tests passed in 0.89 seconds.
These suites overlap and are not additive independent tests. Coverage
includes both arms, wrong/live scheduler records, altered source/context,
mapping/count/payload changes, helper provenance, output inventory, fresh
namespaces, failure retention, native endpoint checks and arithmetic.
Both batch scripts passed `bash -n`; both frozen-worktree CLI imports passed
`--help`; staged whitespace checks passed. No real CPM scoring smoke test
has run because upstream scientific prerequisites remain pending.

Next: extend the corrected SwissTree raw-family count auditor to admitted
CPM outputs, then connect all six parameter arms to the frozen 18-contrast
paired bootstrap. Missing/failed arms must remain explicitly missing,
without multiplicity changes or adaptive retries. These submissions do not
close the parameter-robustness requirement by themselves.
