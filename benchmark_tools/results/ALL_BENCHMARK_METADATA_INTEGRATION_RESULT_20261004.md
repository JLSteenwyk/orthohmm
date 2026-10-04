# Integrated Historical Benchmark Metadata

The [current integrated table](all_benchmark_metadata_integrated_20261004_v3/register.md),
[structured register](all_benchmark_metadata_integrated_20261004_v3/register.json)
and [resource table](all_benchmark_metadata_integrated_20261004_v3/resources.tsv)
combine the existing all-tool register and the separately verified Three
Kingdoms historical metadata supplement. No earlier audit, native inference,
timing or scoring is repeated. Both older inputs remain unchanged.

## What Is Delivered

Exactly eight methods appear on each of OrthoBench, corrected QfO and Three
Kingdoms: 24 method/dataset rows. Every pre-existing row field survives exactly,
including scores, output/pair semantics, declared versions, commands, prediction
and input records, resource scopes, limitations and missing-evidence notes.
The new layer adds supplemental fields only for the seven historical Three
Kingdoms methods. The newer matched SonicParanoid row is untouched; its replaced
historical run is not joined to current scores.

The added fields include both OrthoHMM metrics/recorded argv, six GNU-time
commands and intervals, launch-source versus post-run harness chronology,
inherited manifests, dirty-worktree flags, empty logs and ordered recovery
events. Supplement prediction and version declarations must match the selected
row before integration. Recorded commands or version text do not become new
execution or runtime attestations.

The resource TSV preserves 29 earlier table entries, including unavailable
placeholders, and adds eight historical intervals: 37 entries total, not 37
independent runs. There are 34 numeric wall observations and three explicit
unavailable entries. Two unavailable rows are corrected QfO OrthoHMM's missing
full-inference resources; the other is OrthoBench ProteinOrtho's unestablished
per-invocation wall interval. Unknown values are `NA`, never zero.

Native intervals and integer wrapper observations remain separate. Checkpoint
conversion does not become separately timed sequence-only inference. OrthoMCL
recovered mode4 and BPO conversion do not replace full search-to-clusters cost.
FastOMA retains supplied-tree dependence and driver/container CPU/RSS caveats.
The original register's `full_inference` fields are not rewritten or promoted
by these retained-directory associations.

## Actual Verification

Commit/push the tested reporter at `5c5c0192` before the final collection.
The final v3 invocation exits zero in the installed Python 3.12.3 analysis
environment. Existing imports require Biopython; a minimal `-I -S` invocation
failed before output. No new dependency or environment was installed.

The initial 72-test suite checks the second preliminary export. The source-
frozen preparation suite passes 71 cases and skips one not-yet-created final
export. After v3 collection, all 72 cases pass without skips, including real
metadata/source identities, every selected row, all table fields and scope
refusals. Source and source-pin checks are finalized before the selected run.

Independent standard-library readback, without importing the integration
helper, proves the Cartesian 24-row inventory, all 72 individual score/secondary-
mean positions, unchanged original fields, exact seven supplement joins,
unchanged matched Sonic choice, source hash, all 37 resource entries and all
481 TSV fields. All eight supplemental measurement dictionaries agree with
their retained sources. No resource entry is labelled an independent repeat.

The [execution receipt](benchmark_metadata_integration_execution_20261004.json)
pins the final results, sources/tests, three JUnit files and development evidence.
The first successful preliminary export omitted the two unavailable QfO cells
and had 35 entries; it remains unselected. The second added the placeholders
and had 37 entries, but predates the final refusal/source-binding checks. Both
preliminary source versions are preserved in local byte-verified archives,
not substituted for the final source. No native failure or scientific retry
is inferred from these reporting-development steps.

## Reproduction And Limits

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  benchmark_tools/integrate_benchmark_metadata.py --repo . \
  --output /tmp/orthohmm-integrated-benchmark-metadata
```

Use an absent output directory. The integration directly checks three fixed
metadata/score files, not every transitive raw input, prediction, executable
or historical consumption event. Raw identities/manifests remain explicitly
inherited. Version declarations remain declarations. Scope differences and
unknown, potentially method-dependent shared-host contention forbid isolated
efficiency rankings. This integration does not supply missing full-inference
costs or retroactively standardize historical tool allocations.

The separate matched-resource Threadripper panel and original factorial stage-
cost table remain separate, rather than replacing costs of these score runs.
QfO similarity metrics and its secondary mean are not universal F1. Three
Kingdoms is supplementary BUSCO co-membership, including within-species pairs.

This advances consolidation requirements 1.2/1.4; complete transitive provenance,
historical consumption, unresolved uncertainty, development-family inventory,
biological-stratum coverage, per-configuration full costs and distribution
requirements remain open. The frozen method, scores, main PDF and rc3 archive
are unchanged. The full publication goal is active and incomplete.
