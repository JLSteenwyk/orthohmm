# Corrected SonicParanoid Command Preparation

Prepared, not submitted. The native command preserves the original
SonicParanoid 2.0.9 default invocation: fresh input and output directories,
32 threads, no update/overwrite/reuse option, and no scientific retuning.
The original native log records DIAMOND very-sensitive mode, minimum
bitscore 40, in-paralog length difference 0.75, merging threshold 0.75,
MCL inflation 1.50 and graph-only mode false. Those required settings are
recorded for later native-log admission.

Command manifest SHA-256:
`188efb1f18e8f0bd90967670d18015ad500a9c47480333b0b380dff7c10c4c59`.
Read-only runtime snapshot SHA-256:
`b4ed57defaba417f6627b671df958a2d533f4438bb8ff35a0b95cb820fd61bc9`.
All 78 corrected input files and frozen comparison/primary manifests were
verified. Ten focused tests passed; real manifest preparation succeeded.

## Resolved Dependencies

The probe used the installed resolver without invoking an installer or
inference. All five dependency version probes exited zero:

| Dependency | Observed version |
| --- | --- |
| DIAMOND | 2.1.9 |
| BLASTP and makeblastdb | 2.15.0+ |
| MCL | 14-137 |
| MMseqs | 45111b641859ed0ddd875b94d6fd1aef1a675b7e |

These are SonicParanoid's bundled executables, including MCL under the
non-conda resolver mode. The complete 96-record package inventory and
all five executable identities match the prior biological-application
runtime snapshot. This is current identity evidence, not proof of the
executables used in the historical April QfO run. The default DIAMOND
label does not exclude downstream MMseqs profile searches; the original
native log explicitly records that stage.

The existing full runtime inventories were independently reverified:
132,274 records in `biological_wgd_runtime_trees_v1.json` (SHA-256
`70eda6198d28d6c36c697d2862911a0853d9482c9e15f504f17574ca73f99071`)
and 14,600 in `biological_wgd_system_trees_v1.json` (SHA-256
`5a19ac94e14aff479c7bd7a09a757a8471056d23a9b51fbbc8c51a2d8aee0e27`).
Both matched their exact current trees. They can be bound into prospective
execution and rechecked; this does not authorize inference by itself or
establish a hermetic operating-system snapshot.

## Remaining Gates

No scientific execution is yet authorized. Bind full interpreter/system
inventories, explicit environment and isolated bytecode lookup to a pinned
fresh-copy launcher; recheck effective resolution before/after inference.
Native pair outputs must come from species-to-species ortholog tables,
not clique expansion of global groups. Validate all species-pair tables,
gene ownership and completion before conversion/scoring. Freeze corrected
participant identity and scorer commands before evaluation. This is
shared-host accuracy work, not a dedicated timing run.
