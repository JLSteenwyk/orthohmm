# Fixed-Length Stress Panel: Interpretation

**Superseding runtime warning:** The
[runtime audit](PROFILE_RUNTIME_FAILURE_20260916.md) discovered missing
profile-alignment binaries in the frozen OrthoHMM checkout. OrthoHMM scores
and tree-failure interpretations below require corrected-runtime reruns;
they do not establish behavior of the intended high-sensitivity method.

Source: `simulation_fixed_length_results_20260916.json`, assembled with pinned
scorer a58fd3a after all 70 method-array tasks were terminal. The readable
table is `SIMULATION_FIXED_LENGTH_RESULTS_20260916.md`. The 280 records include
explicit failed outcomes; no failure has an imputed accuracy score.

## Completion and Failure

- OrthoHMM high sensitivity: 70/70 datasets completed and admitted.
- OrthoHMM satellite_v2: 64/70 completed and admitted; six execution failures.
- Full OrthoFinder: 0/70 admitted. Forty-eight runs lack native completion;
  22 have a completion marker but nonfinite graph weights. All returned zero
  at the process level, illustrating why process exit alone was insufficient.
- All 70 OrthoFinder sequence checkpoints are also excluded by the parent's
  native-completion/finite-weight admission rule. Invalid preprocessing is
  not repaired by converting an existing checkpoint file.

The OrthoFinder numerical failure is independently reproduced in
`SIMULATION_NATIVE_FAILURE_AUDIT_20260916.md`. Consequently, all 14 planned
OrthoHMM-versus-OrthoFinder contrasts have no complete paired seeds. There
are no comparative F1 estimates or confidence intervals for those contrasts.
This panel cannot support a claim of general superiority over OrthoFinder.

## OrthoHMM Failure Mode

The six satellite_v2 failures occur at seeds 20261003, 20261006 and 20261007
in both divergent conditions. All six native `benchmark.log` files end with
`PhylogenyConfigurationError`: species-tree inference lacks connected
multi-species single-copy families covering every taxon. The uncovered taxon
is n2 for seeds 20261003/20261006, and n3 plus n4 for seed 20261007.

These are end-to-end failures, not missing-reference exclusions. The frozen
configuration is not changed to make them disappear. The seven successful
seeds in each divergent condition form a selected subset; their satellite_v2
means must not be compared to the ten-seed high-sensitivity means as though
they were a paired improvement estimate. Trace candidate recovery and tree
family selection in subsequent development, with fresh confirmation required
for any changed method.

## Descriptive Results Only

For baseline, all ten seeds succeed for both OrthoHMM configurations: mean
F1 is 96.29% for high sensitivity and 99.68% for satellite_v2. Their precision
and recall trade-offs, all seven conditions, and success denominators are
reported in the generated table. These are descriptive stress-panel results,
not independent curated validation or a matched HMM-contribution ablation.

The heterogeneous-family-length panel was separately frozen before these
accuracy outcomes were assembled. Its generation is complete; all 20 paired
histories match. Exported lengths were rechecked for 33,618 sequence instances
across 40 native exports and all matched their assigned family lengths.
Method inference is ongoing and will be analyzed separately with its own
prespecified seeds and bootstrap settings.
