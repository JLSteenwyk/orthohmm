# Corrected Fixed-Length Simulation Results

[Results table](SIMULATION_FIXED_NATIVE_RESULTS_20260916.md) is generated from
`simulation_fixed_native_results_20260916.json`. Pinned assembler b66225d
verified corrected OrthoHMM array 21143 and original OrthoFinder array 20957
against their own executors, manifests, inputs and native outputs. No failures
were assigned zero scores. Original defective-runtime results remain intact.

## Findings

- OrthoHMM high sensitivity completed and passed admission on all 70 datasets.
- Satellite_v2 completed and passed admission on 64/70. The six failures are
  seeds 20261003, 20261006 and 20261007 in both divergent conditions. The
  corrected native logs still report insufficient connected single-copy
  families to represent n2 (seeds 3/6) or n3 and n4 (seed 7) for species-tree
  inference. Restoring the native profile runtime did not resolve this failure.
- Across all 140 corrected OrthoHMM metrics records, profile counts range
  from 62 to 113; every record reports zero added profile edges. Thus this
  panel executes profile construction but provides no positive evidence that
  profile expansion improves its groups. Initial HMM search is still present.
- All 280 outcome statuses and score dictionaries match the original snapshot
  exactly. This does not retroactively validate the defective runtime; it
  establishes stability of these outcomes under the documented runtime repair.
- No OrthoFinder output passes the original native-completion/finite-graph
  gates. The previously diagnosed constant-length normalization issue remains
  separate from OrthoHMM's runtime defect. No new OrthoFinder run was needed
  or silently substituted. Its sequence checkpoint retains the parent gate.
- All 14 planned OrthoHMM-versus-OrthoFinder contrasts have zero successful
  paired seeds. Differences and intervals are unavailable, not zero or wins.

The panel is a simplified fixed-length stress test, not a general comparative
accuracy result. Means for satellite_v2's divergent conditions are conditional
on successful runs; do not compare them to ten-seed means as paired effects.
Shared-machine runtime logs do not establish matched performance advantages.
Variable-length results must be reported separately; curated validation,
ablations, robustness, application and publication requirements remain open.
