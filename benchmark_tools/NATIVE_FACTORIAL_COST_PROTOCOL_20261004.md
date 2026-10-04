# Remaining Native Factorial Cost Protocol

## Scientific Scope

The [eight-cell adapter diagnostic](results/NATIVE_FACTORIAL_ADAPTER_RESULT_20261004.md)
provides an executable benchmark-only adapter for the frozen method at
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`. P disables only downstream
profile expansion, not initial HMM search, k=4/cap100 sensitivity, multipass
RBNH/singleton grouping, Leiden CPM0.1/seed4 or default refinement. C uses
the existing satellite_v2 policy with unchanged parameters. R selects
fresh MAFFT/FastTree species-tree inference and reconciliation with
species-overlap rooting, positive-paralogy pairs and minimum-variance
species-tree rooting. BLOSUM62 and E1e-4 stay fixed.

For C-on/R-off only, the adapter removes the reconciliation condition from
the existing candidate-expansion gate. This is an explicit experimental
execution path, not a production default or a new accuracy method. Candidate
groups are the R-off output; R-on outputs are root HOGs for OrthoBench and
native resolved cross-species pairs for QfO. Do not treat candidate edges,
clusters and resolved pairs as interchangeable predictions.

This protocol is frozen before any new full-dataset timing under this
adapter. It is not an execution authorization or a passed launch gate.
The bounded diagnostic CLI explicitly refuses full benchmark-sized input.
A separately validated full-run executor and live resource handoff are
still required. Existing successful diagnostics, accounting calibration
and the completed 27-attempt scaling panel must not be repeated solely
because this adapter exists.

## Remaining Identities

Reuse the [two native OrthoBench configuration associations](results/FACTORIAL_NATIVE_RESOURCE_LINKAGE_RESULT_20261004.md)
and [one historical native QfO association](results/QFO_NATIVE_COST_RESULT_20261004.md).
Check any additional retained measurement's exact bindings before releasing
a new identity; do not assume this combined inventory exhausts every
historical output. The following prospective order covers the thirteen
configurations still unassociated in those reports:

| Index | Dataset | Cell | Repeat |
| ---: | --- | --- | ---: |
| 0 | OrthoBench | p0_c0_r0 | 0 |
| 1 | OrthoBench | p0_c0_r1 | 0 |
| 2 | OrthoBench | p0_c1_r0 | 0 |
| 3 | OrthoBench | p0_c1_r1 | 0 |
| 4 | OrthoBench | p1_c0_r1 | 0 |
| 5 | OrthoBench | p1_c1_r0 | 0 |
| 6 | Corrected QfO | p0_c0_r0 | 0 |
| 7 | Corrected QfO | p0_c0_r1 | 0 |
| 8 | Corrected QfO | p0_c1_r0 | 0 |
| 9 | Corrected QfO | p0_c1_r1 | 0 |
| 10 | Corrected QfO | p1_c0_r1 | 0 |
| 11 | Corrected QfO | p1_c1_r0 | 0 |
| 12 | Corrected QfO | p1_c1_r1 | 0 |

One new attempt per missing configuration supplies descriptive costs, not
stable per-cell variance or a fully replicated causal ablation timing
panel. The existing two OrthoBench configurations retain all three repeats
each; historical QfO retains its one observation. Do not rerun them to seek
faster times or silently count a cached stage as a full pipeline.

The twelve OrthoBench FASTA records are pinned by
`orthobench_factorial_prepared_20260916.json` (SHA256
`5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382`).
The 78 corrected QfO FASTAs are pinned by
`benchmarks/results/qfo_corrected_factorial_v1/manifest.json` (SHA256
`d8385c50426e690afd6d32f3c5302e678de6977c9841451d013201b0f75b564a`).
Rehash the full input manifest before/after each attempt, freeze and verify
native file enumeration/ID ownership, and do not replace these bytes with
historical original-release or sanitized comparator inputs.

## Execution And Accounting

Use the local Threadripper only. Freeze 32 distinct native CPUs0-31, reserve
the corresponding64 Slurm SMT slots, request/enforce128GiB RAM and use
the existing private x86 Python3.10 scientific deployment. Preserve the
single-threaded numerical-library settings, fixed hash seed, fresh absent
bytecode prefix and explicitly bound frozen package/compiled helpers.
Run sequentially without node exclusivity or interference with unrelated
work. Native search uses8 workers times4 threads. The prospective native
timeout is85,800s with the existing26-hour scheduler envelope; report timeout
as a retained attempt, not a reason to select another outcome.

Reuse valid retained runtime/accounting/observer calibration bindings. The
new executor must separately bind its own source, adapter, complete inputs,
native command/factors, private-runtime lookup and measurement workflow.
It must verify current scheduler/cgroup CPU/RAM/affinity limits and actual
safe capacity at the native handoff. Contention alone is accepted; do not
introduce a quiet-window, DGX or observer-isolation gate. Do not stop,
suspend, renice, re-affinitize or alter unrelated jobs or services.

Use the existing lifetime native-step cgroup memory peak and native subtree
CPU-stat bracket with a monotonic launch-to-exit interval; preserve their
exact launcher/wrapper scopes. Keep runtime/input hashing and materialization
outside the native interval, and report their costs separately. Fresh full
inference must not load a historical hit/graph/profile/tree checkpoint.
Its own newly produced checkpoint is retained for verification.

Observe whole-run foreign CPU demand, memory/swap pressure and relevant I/O
using the retained periodic observer and scope-attribution policy. Record
the launch inventory, cadence violations and all terminal evidence. Matching
resource limits does not establish isolated performance. Distortion is
unknown and may differ between configurations; no background subtraction,
corrected isolated estimate or definitive speed rank is permitted.

## Admission And Reporting

Retain every attempt, failure and timeout under a unique identity. Never
overwrite output, implicitly resume a failed identity, retry until faster,
select a fastest repeat or fill a failed cell from another execution. A
genuinely necessary replacement requires a separately prospective amendment
and keeps the original attempt visible.

Before release, verify the preceding identity's terminal state from Slurm
and actual measurement/native handles. After completion, independently
review source/runtime, output universe/factor/stage/checkpoint semantics,
resource accounting, environment and scientific conversion bindings.
Compare resulting whole partitions against the original factorial outputs;
retain differences and distinguish configuration-level association from
exact output reproduction. Do not select only partition-matching outcomes.
Accuracy equivalence is not implied by matching a gene universe or counts.

Report each observed full native interval and every unavailable cost.
The original sixteen cached factorial full-cost fields remain null; later
native observations are separate, not retroactive measurements of those
executions. Keep old GNU-time/process-tree RSS, current lifetime cgroup peak
and original/private deployments separate. Existing mixed-deployment costs
are descriptive only and do not support matched causal component overhead.
New native output changes must not inherit old accuracy claims without
appropriately scoped output/conversion/scoring validation.

These new costs do not close broader independent-validation, QfO uncertainty,
prespecified biological strata, all-method restoration or public archival
requirements. Publication readiness remains unproven.
