# Native OrthoBench Factorial Scoring

## Scope

The separate [scorer](score_native_orthobench_attempt.py) evaluates completed
full-native OrthoBench cells against all70 frozen RefOGs. It does not change
any source bound to the currently running22428, rerun inference or release a
subsequent identity. This completes the scoring workflow, not the live cell
or wider publication goal. Corrected QfO still requires explicit pair
conversion and its separately validated native six-endpoint assessment; a
RefOG statistic is never substituted for those endpoints.

## Admission

Require the exact request/plan/index/cell/dataset/repeat and a successful
terminal join of runtime, resources, environment and semantic outputs. Query
Slurm again before scoring; a live controller must also retain the request
digest in its comment. Check all current bound helper/input pins and the
output validator's checked artifacts. Reject incomplete, live, failed,
timed-out or differently bound attempts before writing a score. Separately
recovered22427 remains scientifically usable via its explicit recovery
receipt; this route never relabels it as a clean native success.

The [actual live refusal](results/native_orthobench_live_refusal_22428.json)
uses the real22428 request/plan and retained startup observation, with a fresh
squeue query confirming RUNNING. The scorer refuses it and creates neither a
score nor destination. This is an admission-boundary check, not an accuracy
result or terminal runtime review.

## Statistic And Output Meaning

R-off cells use final space-separated clusters; R-on cells use the native
root-HOG TSV. Both are **group co-membership** endpoints here, not a claim
that every within-group pair is a resolved ortholog. Actual native pair
orthology is reserved for the corresponding pairwise evaluation workflow.

Apply the existing audited official-formula scorer, all70 frozen reference
families and their low-certainty exclusions. Independently reproduce the
earlier cached cell's F1/precision/recall before comparison. Recompute the
weighted statistic from aggregate weighted TP/FP/FN; do not average family
F1 or silently substitute a metric. Require identical full protein universes
and report complete canonical-partition differences, reference coverage,
per-family sufficient statistics and cached/native score deltas. Changes are
retained rather than treated as automatic failures or parameter-selection
signals. These comparisons establish reproducibility, not a new method effect.

Store an immutable JSON score with request/plan/review, fresh scheduler,
source/reference/prediction provenance and checked artifacts. Scientific
development exposure and lack of independent validation remain explicit.
The score links to observed shared-host resource scopes without substituting
internal stage time, correcting contention or claiming isolated efficiency.
No bootstrap uncertainty or next-job authorization is implied by scoring.

## After Actual Terminal State

First apply `review_native_factorial_attempt.py` to the same22428 request,
using the matching Python3.10 review environment and a new, separate review
directory. Its retained failure is not retried merely for a faster/clean time.
Only if independent terminal review returns native_success, run:

```bash
benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/score_native_orthobench_attempt.py \
  --request benchmarks/work/native_factorial_launch_20261004/request_01_receipt_amended.json \
  --request-sha256 ba72ed12e822465140257b0fe15acc56dbe12fd8ef8a398abae075afe00df02d \
  --terminal-review /absolute/path/to/review.json \
  --terminal-review-sha256 ACTUAL_RETAINED_REVIEW_SHA256 \
  --output-directory /absolute/path/to/new/scoring_directory
```

Use the sanitized scientific review environment as in retained workflows;
never invent the future review hash or reuse a score directory. The output
does not execute QfO, change resource/source plans or release index2.

## Verification

[Final616-test report](results/native_orthobench_scoring_tests_20261004.xml)
passes with zero failures/errors/skips in10.25s. It includes44 new tests:
terminal/identity/resource admission, all8 factor combinations' format
selection, weighted versus averaged F1, low-certainty exclusions, full
universe/partition checks, changed cached points/artifacts, joined successful
scoring, fresh terminal/comment failures, resealed wrong-output rejection and
no-overwrite behavior. Joined scheduler/source tests use explicit doubles and
synthetic inputs; they do not certify an actual completed native cell.
An intermediate2failure/42pass run lacked the input-list field in its joined
synthetic fixture; fix the fixture, not production requirements. Its report
remains in work. No search, phylogeny or diagnostic inference was rerun.
