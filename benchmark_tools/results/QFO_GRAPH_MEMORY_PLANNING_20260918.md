# Sequence-Control Graph Memory Planning

## Finding

The frozen graph code memory-maps numeric checkpoint inputs but subsequently
allocates full-length boolean masks, selected hit copies, integer slot keys
and winner arrays. Therefore checkpoint loading is not evidence that the full
all-hit graph fits a given memory allocation. Do not silently truncate the
all-hit variant, and do not treat the independently prespecified top100
diagnostic as a replacement if all-hit execution fails.

Frozen `orthohmm/accuracy.py` SHA-256:
`1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6`,
from core revision `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`.

## Allocation Calculation

`estimate_rbnh_array_payload.py` audits the checkpoint and verifies these
frozen source bytes. For finite native int32/int32/float64 hit arrays, let:

- H = total directed hits, E = nonself hits.
- S = genes times (maximum species ID plus one), matching native allocation.
- W = tied winning hit rows and K = occupied query/target-species slots.

At the source line immediately after allocation of `best_targets`, named
array payload is exactly `H + 24*E + 12*S + 8*W + 24*K` bytes. The tool
bounds this snapshot using `1 <= K <= W <= E` and `K <= S` when E is positive.
When no eligible hit exists, native code returns before this snapshot.
Input array logical size is reported separately as `16*H + 4*genes` bytes.
The model requires 64-bit NumPy indices. It does not measure or predict RSS.

Excluded costs include sorting temporaries, allocator overhead, gene strings,
input page residency, libraries, later edge deduplication, graph-library
construction, clustering, singleton assignment and refinement. In particular,
the snapshot upper bound is NOT an upper bound on total peak RAM. There is no
automatic graph-feasibility admission or recommended memory allocation.

## Completed OrthoBench Check

Both already-retained checkpoints were rehashed and numerically audited;
no graph inference was repeated. File identities were checked again after
the calculation. These observations are not extrapolated to QfO.

| Variant | Directed hits | Snapshot payload bounds, bytes | Input logical bytes, separate |
| --- | ---: | ---: | ---: |
| All hits | 100,099,147 | 2,532,653,619 to 3,403,835,787 | 1,602,591,864 |
| Post-search top100 | 48,991,663 | 1,254,966,519 to 1,717,288,815 | 784,872,120 |

Machine-readable reports: `ob_all_hits_rbnh_array_payload_20260918.json` and
`ob_top100_rbnh_array_payload_20260918.json`. Each binds its checkpoint,
numeric audit, individual input files, frozen core and estimator source.
27 focused tests pass, including eight randomized tied-score cases that
trace actual native named-array sizes at the specified line, early return,
large integer counts, malformed dimensions, wrong core bytes and a checkpoint
change after numeric audit.

## Corrected QfO Next Step

### Queued Admission-Gated Calculation

Job21798 is queued after successful numeric admission21792; the scheduler
confirms PENDING(Dependency). The frozen executor is
`928eba0430f2348647d30819d5a6fa23873cb19d` at
`benchmarks/work/publication_qfo_graph_payload_v1`. It requests2CPUs,64GiB,
4hours onbizon with no requeue. No graph inference is launched.

`plan_qfo_sequence_graph_memory.py` requires completed numeric-admission and
conversion accounting, exact frozen converter/admitter revisions, successful
source-equivalence status, both prespecified variant identities and the
recorded evidence file hashes. It runs the existing estimator against each
admitted checkpoint manifest hash and compares genes/hits/self-hit counts
with the numeric audit. It rechecks evidence before writing the report.
The estimator verifies the frozen native graph source and checkpoint contents.

Future output: `benchmarks/work/qfo_graph_payload_20260918.json`. An existing
output is refused. The report records both variant estimates and explicitly
sets graph feasibility, graph launch, accuracy evaluation and publication
readiness to false. A resource decision still requires review of actual
counts and costs beyond the modeled snapshot; top100 cannot replace all-hit
failure. No completed QfO payload calculation is claimed at submission.

42 focused tests pass, including8 new tests for both variants/exact hashes,
missing variant, path/count disagreements, false feasibility admission and
failed scheduler rejection before output creation. Batch shell syntax passes.

### QfO-Specific Graph Plan Builder

`prepare_qfo_sequence_graph.py` prepares both graph commands after the payload
job completes and its results are reviewed. It requires an explicitly supplied
payload-report SHA-256, successful payload/admission accounting, the frozen
planner revision, unchanged recorded source evidence, and separate integer
memory allocations for all-hit and top100. It refuses an allocation that
cannot hold even the modeled lower snapshot bound; passing this necessary
condition is not proof that the complete graph fits.

Each variant binds the admitted checkpoint manifest and gene-order file,
984,137 genes/78 species,32CPUs,BLOSUM62,CPM0.1,Leiden seed4, profile-off
replay, and two expected clustering calls (`initial`, `multipass`). The
runtime is checked against the frozen core/launcher. Candidate expansion
and reconciliation remain off. It writes only a command plan, not predictions.

The old `run_sequence_graph_control.py` entry point must not be repurposed
unchanged: its admission paths, scheduler IDs and count validation are specific
to OrthoBench, and it predates the checked QfO clustering constructor path.
Only its pure `graph_command` helper is reused. The new plan explicitly
requires a separately frozen checked-clustering executor and prohibits direct
execution of its native command. The checked executor and post-run admission
are now implemented below; actual execution/admission remain pending.

44 focused tests pass (16 new plan-validation/command tests,8 payload-wrapper
tests and20 existing estimator tests); the CLI help smoke check passes.
No actual command plan, allocation decision, graph submission or accuracy
result has been produced while the prerequisite payload report is pending.

### Checked Sequence Executor

`run_qfo_sequence_graph.py` uses the frozen profile-off replay with the
existing checked Python-pair constructor and native optimizer observations.
`sequence_graph_evidence.py` validates the exact plan, variant, payload
estimate, checkpoint and gene-order records. Sequence provenance is a distinct,
mutually exclusive worker mode, never fabricated corrected-HMM admission.
Only `initial` and `multipass` clustering stages are accepted. Historical and
corrected-HMM worker modes remain supported.

The driver requires a scheduled32CPU allocation with memory matching the
reviewed variant plan, unchanged frozen runtime, and a fresh output directory.
It retains payloads, native observations, GNU-time/log output and failure
reports. Both checked calls, profile-off metrics, planned command/paths,
complete gene coverage and postflight runtime/source identities are required
before reporting completion pending admission. The batch entry point
`qfo_sequence_graph_batch_20260918.sh` requires explicit `--mem` at submission.
No graph job is submitted yet and no actual QfO result is claimed.

98 focused tests pass, including real small igraph/Leiden fixtures under
historical, corrected-HMM and sequence provenance, corrupted boundaries,
constructor inputs, modules, partition coverage, retained outputs, wrong stage
counts and mixed/incomplete provenance. A fresh-process smoke check against
the actual frozen launcher confirms correct scientific-module/helper paths.
Full-scale execution and actual post-run admission remain unverified.

### Post-Run Admission

`admit_qfo_sequence_graph.py` requires a terminal successful32CPU job with
the planned memory allocation before reading run outputs. It binds the plan
and parent result to supplied hashes, rechecks numeric checkpoint contents,
scientific settings, source files, runtime observations and exact variant
identity. Both stage partitions must cover the complete corrected universe.
It repeats validation of retained native constructor/optimizer observations
and verifies the multipass output came from the corresponding checked stage.
Altered commands, stage order, counts, file identities or memberships fail.

The executor is pinned to69adea4d5a1dc634783b20ff7bbe9fbe1c6464db in
`benchmarks/work/publication_qfo_sequence_graph_v1`. Generate the eventual
command plan with that frozen executor's preparer so its source record agrees.
The admission script is frozen atf1e21b09c28f270dc3ef2243bdcad86f212b58a0 in
`benchmarks/work/publication_qfo_sequence_graph_admission_v1`. All-hit and
top100 have separate admission reports; neither substitutes for a failed arm.

154 focused tests pass, including25 new parent/partition/gating tests and
expanded stage auditing for both sequence variants. Actual small native
igraph/Leiden corruption tests remain included. The CLI smoke check passes.
This is implementation/test evidence, not admission of a QfO run. Preserved
native observations are checked, not an independent historical memory trace
or a new graph inference. Pair conversion, accuracy scoring, paired uncertainty
and dedicated timing remain separate requirements.

### Group-To-Pair Conversion

`prepare_qfo_sequence_pairs.py` accepts only a successful, source-pinned
sequence graph admission for the selected variant and the final
`multipass_refined` partition. It requires the admission job to have completed
under2CPUs/192GiB and requires its source to match the frozen admission
worktree above. The corrected search input plan, complete FASTA inventory,
prediction, QfO mapping and transitive recorded evidence are checked.

The existing orthogroup converter emits cross-species clique pairs with bare
accessions. An independent per-group species-count formula checks the expected
pair count. Unexpected mapping loss, count disagreement, duplicate identifiers,
unknown genes or incomplete input coverage fail. A valid zero-pair prediction
is retained, not discarded based on accuracy. Partial output and failures are
preserved; final pair files appear only after successful checks.

The report labels these as group-derived clique pairs, never native
phylogenetic ortholog predictions. Future outputs are separate directories in
`benchmarks/results/qfo_sequence_pairs_v1/{all_hits,top100}`. The batch entry
point `qfo_sequence_pairs_batch_20260918.sh` requests2CPUs/64GiB/24hours;
it has not been submitted.33 focused tests pass, including actual converter
subprocesses and corruption controls; CLI smoke and batch syntax checks pass.
QfO assessment and score admission are connected below, but have not executed
on actual sequence-control results.

### Endpoint Assessment And Admission

`run_qfo_sequence_assessment.py` validates the terminal2CPU/64GiB conversion,
the frozen converter, variant identity, pair semantics/counts, zero mapping
loss, mapping identity and frozen assessment environment. It uses the existing
six-endpoint command builder, with separate work/result directories for each
variant. Zero-pair conversions are not rejected on score grounds; a resulting
native scoring failure is retained without retry or imputed scores.

The converter is frozen at4f0c30e5cdf287a35c9600886aec0a41bcc0b720 in
`benchmarks/work/publication_qfo_sequence_pairs_v1`. The scoring runner is
frozen at425f0a7f5d5dc9e1438ab0a1766c45b13596dc14 in
`benchmarks/work/publication_qfo_sequence_assessment_v1`. Its batch entry
point `qfo_sequence_assessment_batch_20260918.sh` requests8CPUs/64GiB/24hours.
Neither variant is submitted while its graph and conversion results are absent.

`admit_qfo_sequence_assessment.py` independently requires terminal scoring
success with the requested resources, reconstructs expected execution records
using helpers identical to the frozen runner, and checks preflight consistency,
complete output inventory, native task trace, and all six endpoint files using
the existing native-assessment validator. Process success alone never admits
accuracy. Existing outputs cannot be overwritten. The endpoint mean remains a
project-defined secondary summary, not official QfO F1, matched sensitivity,
paired uncertainty or independent biological validation.

56 runner-related tests and74 admission/native-validator tests pass (overlapping
coverage, not130 distinct tests). CLI and batch syntax checks pass. Actual
QfO sequence-control graph execution, conversion, six scores and paired
uncertainty remain pending; implementation is not evidence of those results.

### Earlier Search Progress

As of the recorded live check, DIAMOND job21789 remained running; its
incremental execution log recorded56 completed target searches out of78,
27,688,723,204 bytes of completed hit tables, and active zero-based target56.
These are progress observations, not a full-panel success claim or hit count.
Do not estimate the remaining workload by assuming equal-sized targets.

After independent numeric admission21792, apply the estimator separately to
the admitted all-hit and top100 checkpoints using their exact manifest hashes.
Review the actual numeric counts and file sizes together with available
resources and downstream graph costs before freezing the graph launcher.
The corrected HMM/DIAMOND hit-coverage diagnostic21793 remains separately
queued; it does not establish graph memory feasibility or matched sensitivity.
No graph job is submitted by this planning work, no running pipeline is
changed, and no new accuracy result is admitted.

Example after numeric admission, substituting an admitted checkpoint/hash
and a fresh output path:

```bash
python benchmark_tools/estimate_rbnh_array_payload.py --checkpoint CHECKPOINT --sha256 MANIFEST_SHA256 --core benchmarks/work/publication_method_native_v2/orthohmm/accuracy.py --output NEW_REPORT.json
```
