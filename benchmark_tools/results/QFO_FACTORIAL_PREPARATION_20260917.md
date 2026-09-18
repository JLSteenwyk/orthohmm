# QfO Factorial Candidate Preparation

Preparation job21670 completed successfully (exit0:0, elapsed00:03:04) on
bizon. Its pinned executor was bd5229db0b9e37396a8dd54f463c2ccd743c2527.
The retained preparation manifest is `qfo_factorial_prepared_20260917.json`,
SHA256 `706b07c91e9a130dae229837641a7ad62d7d36a09679e5e0daa0959e182b7d64`.
Original artifacts remain at `benchmarks/results/qfo_factorial_v1`.

All four candidate arms preserve the validated 976,504-gene universe and
78-proteome ownership. Expansion uses the frozen satellite_v2 engine:

- Profile-off: 393,231 seed families; 40,482 merges; 352,749 candidates.
- Profile-on: 390,980 seed families; 40,073 merges; 350,907 candidates.
- Both expanded arms complete two iterations and retain native seed-membership
  constraints. The unexpanded arms preserve their corresponding seed partition.

These are grouping/preparation counts, not accuracy outcomes. Native core,
compiled runtime and numeric checkpoint checks agree before and after
preparation. The output records incremental preparation costs only, not
controlled end-to-end inference performance. The original failed submission
21669 remains a commit-argument gate failure before Python execution.

## Reconciliation Gate

`run_qfo_factorial_cell.py` accepts only the complete hash-pinned manifest and
a matching terminal successful preparation job. It regenerates all eight
planned commands, checks input/candidate/constraint/source hashes, validates
each seed against recovered admission, and verifies frozen runtime and tool
environment before and after processing. All four cells passed check-only
validation without inference or reference scoring.

The planned executor contains newer development core code and must NOT be the
reconciliation import root. The runner relocates only the launcher path to
`publication_qfo_replay_native_v1`, after proving the launcher and its imported
benchmark helper byte-identical to the prepared sources. Scientific arguments,
candidate paths, destination paths and interpreter are unchanged. The actual
command and both source identities are recorded. The native core/runtime
verification independently binds that launcher to the frozen publication core.

Reconciliation is scheduled sequentially on bizon with32 CPUs/192 GiB, up to
48 hours per cell, no automatic requeue or artifact overwrite. The DGX timing
allocation is untouched. The process wrapper captures GNU time and hashes
outputs; these measurements remain shared-workstation incremental costs.
Successful process return is not native admission: complete RootHOG coverage,
reconciliation manifests, source-family constraints, conversion and native
scoring remain required. No new accuracy scores are available yet.
