# Search Descriptions Corrected Without Numerical Changes

Corrected comments/docstrings in six current search modules against the
[frozen numerical specification](FROZEN_HMM_SCORING_SPECIFICATION_20260930.md).
The frozen revision and all historical source/output identities remain intact.
No benchmark inference, conversion, score, default or threshold was changed.

## Corrections

- Single-sequence profiles use integer substitution rows, uniform insert
  emissions of -1 and additive transitions, not normalized HMMER probabilities.
  The background-frequency argument does not change these emissions.
- Default closing cost is -1, not the previous comment's -2. Transition storage
  is (7,), not (L, 7). Penalties affect the raw-score significance gate before
  length normalization; normalization does not make them irrelevant.
- The local match recurrence includes its zero restart, and the final score
  is bounded below by zero. This is maximum-path scoring, not Forward scoring
  or an orthology posterior. Exact Plan7/phmmer equivalence is not established.
- The static E-value approximation is not demonstrated calibration for the
  actual gap costs and banded recurrence. Stored initial normalized scores
  are not calibrated HMMER bit scores. Matrix constants are unchanged.
- MSA columns at exactly the gap/unknown cutoff are dropped. Consensus uses
  the maximum pseudocount-adjusted probability, not the raw majority residue.

## Validation

The [equivalence receipt](search_description_equivalence_20260930.json) binds
the six before/after file hashes to parent revision
`4fbae61713149ff589fc0b0bb05af814c2b9ee2a`. Its SHA256 is
`63f388a0c652576d44bdfd80286901e05134a1f461435f580c09f2a2720e6395`.
System Python read the parent blobs with `git show`, parsed each source using
`ast.parse`, removed only leading module/class/function docstring statements,
and compared `ast.dump(include_attributes=False)`. All six match exactly.
Comments/source locations are absent from this comparison; documentation
metadata intentionally changes. This is not proof of all runtime/backend
equivalence, empirical significance calibration or native memory safety.

Focused current-code validation completed with exit 0:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -B -m pytest -q \
  tests/unit/test_profile.py tests/unit/test_evalue.py tests/unit/test_matrices.py \
  tests/unit/test_engine.py tests/unit/test_engine_score_routing.py \
  tests/unit/test_msa_native_dependency.py
```

Result: **130 passed in 3.74s**. This test duration is not controlled performance
evidence. `git diff --check` also passes. The unchanged frozen readback was
reused, not rerun. No completed scientific workload or diagnostic was repeated.

## Remaining Scope

Documentation corrections do not admit a new scientific method, establish
calibration, repair high-CPM admission 22155 or prove publication readiness.
Historical source pins remain historical; a future final timing/release
recipe must bind its actual source bytes. Existing manuscript archives retain
their own original contents and have not been regenerated for this change.

Retained [GDB diagnostic 22159](QFO_CPM_BACKTRACE_RESULT_22159.md) already failed
to reproduce the crash. The isolated reader and forced-GC controls likewise
do not identify its cause. No new debugger run or admission retry was launched.
The Threadripper timing panel remains deferred, without a new contention probe,
quiet-window request, DGX access or unrelated process/service action.
