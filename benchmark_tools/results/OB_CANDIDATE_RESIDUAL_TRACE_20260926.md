# Residual Candidate Difference

The [machine-readable trace](ob_candidate_residual_trace_20260926.json)
compares the validated job-22320 0.11.0 arm to retained historical outputs.
Three source reports are hash-pinned; all selected partition and merge-trace
files are checked against their recorded identities before and after reading.
This is a post-hoc localization diagnostic, with no inference or scoring.

All four pre-candidate outputs are **byte-identical**, not merely equal in
group count or invariant to labels:

| Retained stage | Groups | Historical/replay agreement |
|---|---:|---|
| Multipass | 51,181 | Exact bytes |
| Multipass refined | 63,245 | Exact bytes |
| Profile-final | 50,894 | Exact bytes |
| Profile-final refined | 62,885 | Exact bytes |

The last partition is the explicitly verified seed for historical candidate
expansion. Both runs use the same recorded candidate parameters. Candidate
expansion is therefore the first observed partition difference: both have
54,445 groups, but 12 groups per side differ, involving 93 genes. All changed
groups and genes are retained in the JSON, not selected by benchmark labels.

Both traces contain 8,440 merges. Comparing iteration and directed source/
target gene sets gives 8,427 shared events and 13 distinct events per side.
Numeric cluster labels and within-group gene order are ignored for this
semantic comparison. The first event-order difference is at zero-based index
91; it is not established as the first causally consequential event.

Among shared events, support differs for 1,661 events with maximum absolute
difference 2.1316282072803006e-14; margin differs for 1,520 events with maximum
absolute difference 5.684341886080802e-14. These observations motivate a
numerical-ordering hypothesis, but do not test it. Candidate selection uses
support/margin ordering and bounded attachment counts; inferring causality
from small numerical differences alone would be premature.

The next bounded intervention should hold the exact seed partition, gene
indexing, species mapping, frozen scientific source and runtime constant,
and independently vary historical/fresh normalized score values and directed
hit order. Self hits must be handled explicitly so they do not become an
unreported extra factor. Preserve every arm; do not round scores, change
tie-breaking, tune thresholds, or rerun full search/phylogeny at this stage.
Only after that comparison should a reproducibility fix or further native
validation be considered. Final F1 causality remains unresolved.

Fifty-three focused tests pass for the trace, dependency reader and strict
candidate auditor. Tests distinguish reordered events from changed event
membership, ignore irrelevant labels, and reject malformed/duplicate merges
and NaN values. Historical scores/defaults remain unchanged.

```bash
python -m benchmark_tools.trace_ob_candidate_residual --repo . \
  --output /tmp/ob-candidate-residual.json
```
