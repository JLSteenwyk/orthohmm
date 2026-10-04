# Selected QfO OrthoHMM Stage Linkage

## Scope

Join the two selected corrected QfO OrthoHMM score rows to their exact
conversion, candidate admission, candidate preparation, admitted checked
replay and, for `p1_c1_r1`, native pair/reconciliation metadata. This is
reporting from retained evidence, not another native run, conversion, scoring
run or complete transitive provenance audit.

The selected high-sensitivity cell is `p1_c0_r0`; its pairs are group-derived.
The selected satellite phylogenetic cell is `p1_c1_r1`; its pairs are native
reconciliation predictions. Both use the same admitted strict-profile seed
and checked cached replay. Preserve that shared identity rather than treating
the two replay associations as independent repetitions.

The fixed inputs are the previous integrated 24-cell register, the corrected
eight-method QfO manifest and the retained 16-cell factorial stage-cost table.
Their checksums are frozen in `link_qfo_orthohmm_stages.py`. Directly check the
selected score admissions and metadata chains. Reuse prior raw-input/output,
checkpoint, runtime and source admissions; do not repeat their expensive
inventories or assert new raw scientific/accounting validation.

## Execution

Prepare and push the tested source before exporting the selected result. Use
the existing scientific analysis environment because the existing provenance
reader imports Biopython. No package installation or Slurm job is needed.

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  benchmark_tools/link_qfo_orthohmm_stages.py --repo . \
  --output benchmark_tools/results/qfo_orthohmm_stage_metadata_20261004
```

The command refuses an existing destination, including a broken symlink.
Retain every original register field and every earlier resource entry. Add
one `qfo_stage_provenance` field to the two selected QfO rows only. Preserve
full-pipeline unavailable placeholders and all score positions. The combined
resource TSV should have 42 entries: 37 earlier entries plus five stage
associations, with 39 numeric walls and three unavailable entries. These are
not 42 independent timed runs or five new measurements.

## Validation And Limits

Portable synthetic tests cover score, cell, prediction, admission, replay
command/checkpoint/seed, input, cost, constraint and native-output mismatches;
invalid numeric measurements; unknown-versus-zero semantics; and refusal of
existing outputs. The actual output readback checks every original field and
all TSV fields. The real source-frozen export must also be independently
read back without importing the reporter: check the 16 direct metadata pins,
both score/conversion chains, distinct interval scopes and resource values.

Cached replay excludes the initial HMM search. Candidate preparation excludes
cache loading and seed construction; CPU and memory were not measured for
those intervals. Reconciliation excludes earlier stages, conversion and
scoring. Internal replay wall/timings and the checked worker GNU-time interval
have different scopes and must not be interchanged. GNU-time process RSS
and sampled summed tree RSS remain distinct. Do not sum stage walls or peaks,
substitute scheduler elapsed time, or impute a whole-pipeline observation.

Any upstream complete native run with partition equivalence remains a
separate observation unless a separately tested linkage establishes its
appropriate scope. No isolated timing ranking, causal overhead, new accuracy,
universal generalization or publication-readiness claim is authorized here.
Contention is unknown and potentially method dependent. No DGX/quiet-window
gate or change to unrelated jobs or services is introduced.
