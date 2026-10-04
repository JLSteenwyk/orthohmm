# Native Factorial Adapter Validated On Bounded Fixtures

## Execution Scope

The benchmark-only [adapter](../native_factorial_adapter.py) supplies the
missing production-independent P-off and C-on/R-off experimental paths
without writing or rebinding the frozen production package. Its source
requires exact pipeline SHA256
`2afb89b9dc683e64e58208188f720e07701d4760ff09d1ac3a53c7a8075b84bb`.
P-off replaces only the immutable accuracy profile's `profile_expansion`
field. The sensitive initial HMM search, k4/cap100, multipass graph and
seed4 stay intact. For C-on/R-off only, AST transformation removes the
reconciliation condition from the existing candidate-expansion gate.
All other statements in `_execute` retain identical AST structure.
Copied function namespaces avoid rebinding production globals.

Prepare/push adapter, bounded driver, probe and initial32 tests at
`56fcc355` before actual native diagnostic execution. The
[probe](native_factorial_adapter_diagnostic_20261004/probe.json) records
**eight successful fresh native diagnostics**, one per P/C/R combination.
The deterministic fixture contains four FASTAs,26 proteins, six designed
sequence families and two extra copies; it is a mechanical fixture, not
a validated evolutionary accuracy simulation.

Use the retained private Python3.10 scientific deployment, the byte-checked
frozen source inventory and pinned MAFFT/FastTree entrypoints. Two inherited
CPU affinity IDs0/1 and2 native threads bound this small diagnostic, not
the production32-core timing allocation. About753GiB RAM is available at
launch; no Slurm job is live. No unrelated job/service is changed.

| Cell | Profiles built | Candidate merges | Reconciled families | Inferred-tree families |
| --- | ---: | ---: | ---: | ---: |
| p0_c0_r0 | off | off | off | off |
| p0_c0_r1 | off | off | 1 | 5 |
| p0_c1_r0 | off | 0 | off | off |
| p0_c1_r1 | off | 0 | 1 | 5 |
| p1_c0_r0 | 6 | off | off | off |
| p1_c0_r1 | 6 | off | 1 | 5 |
| p1_c1_r0 | 6 | 0 | off | off |
| p1_c1_r1 | 6 | 0 | 1 | 5 |

Every cell's fresh initial search produces116 significant hits. All eight
numeric checkpoint manifests and their actual arrays are byte-identical.
P-on builds six profiles and records26 profile hits; R-on genuinely infers
a species tree using five families and reconciles one duplicate-containing
family, recording two duplications. Reconciliation/species-tree checkpoint
reuse counts are zero. Stage inventories exactly follow the declared factors.
All eight final group partitions contain the same26 genes in six groups.
This fixture does not show an accuracy advantage or a profile-induced
membership change: its profile-derived added-edge count is zero.

## Candidate-Gate And Independent Checks

The native sequence fixture has zero satellite merges, so it alone cannot
validate nonzero C-on behavior. Add a separate19-gene control using the
actual frozen candidate block/engine: nine singleton satellites, a ten-gene
anchor and180 directed hits. All eight factor combinations are exercised.
Every C-on block performs eight capped merges in two rounds and leaves
exactly one unattached satellite, **IDg8**, while C-off leaves the seed
partition unchanged. This includes both C-on/R-off paths; it is a controlled
block test, not a second full native benchmark or historical causal replay.

Initial32 adapter cases pass. After actual execution,40 adapter/gate cases
pass. Independent raw-artifact tests plus existing pipeline/accuracy tests
pass58 cases; adding full copied-inventory readback gives **59 passing**,
no failures/skips, in1.57s. Tests check exact allowed AST changes, preserved
production resolver/profile state, separate namespaces, changed-recipe
refusal, all copied raw checkpoint hashes/arrays, native stage/factor data,
input/group universes and actual inferred/reconciled tree outputs.
They do not invoke another native diagnostic or a timing run.

Retain the complete small diagnostic tree: **222 files,363,122 bytes**,
with fixture FASTAs, worker logs, metrics, numeric checkpoints, alignments,
trees and native outputs. The [inventory](native_factorial_adapter_artifacts_20261004.json)
has288 file/directory records and no external symlinks. The unchanged tree
fingerprinter verifies all records; independent tests rehash all222 file
contents and check exact inventory completeness. Absolute paths in retained
JSON describe original local execution; tests resolve copied fixture content
by relative suffix. This is not a hermetic runtime/package archive.
Raw progress logs and trailing empty TSV fields retain their native bytes;
Git whitespace warnings on those artifacts are not normalized away. Source
and documentation pass the separately scoped whitespace check.

The GNU-time [companion](native_factorial_adapter_diagnostic_20261004/companion.time.txt)
records50.47s for the entire eight-cell diagnostic harness and maximum
process RSS218,828KiB. This is not a publication cost or a measurement of
any full benchmark/configuration. It is not added to the scaling panel.
All source, retained probe/inventory and local JUnit bindings are recorded
in the [execution receipt](native_factorial_adapter_execution_20261004.json).

## Next Execution

The [prospective remaining-cost protocol](../NATIVE_FACTORIAL_COST_PROTOCOL_20261004.md)
freezes the13 missing configuration identities/order, exact original input
manifests, factors, shared-host resource/reporting scope and failure rules.
The bounded diagnostic CLI does not authorize full-dataset execution.
Next implement/validate its separately bound full-run executor and actual
native resource handoff, reusing valid existing calibration/runtime receipts.
Then release only genuinely remaining identities sequentially after safe
capacity checks, without a quiet-window/DGX prerequisite or fastest retries.

Full-dataset factor/output equivalence and costs remain unmeasured. No
frozen default, official score, completed scaling measurement, previous
manuscript/PDF/rc4 snapshot, unrelated sample or other workload is changed.
This later diagnostic addendum is not already included in rc4. Broader
scientific validation/uncertainty/strata/restoration/distribution requirements
remain, and the full goal is active with publication readiness unproven.
