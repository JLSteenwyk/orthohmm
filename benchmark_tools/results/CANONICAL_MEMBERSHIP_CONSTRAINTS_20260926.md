# Canonical Membership Constraints

The [constraint audit](canonical_membership_constraint_audit_20260926.json)
compares all five independently admitted canonical candidate arms with the
historical candidate merge trace. Historical preparation, canonical readback,
individual trace files and the reviewed membership-consumer source are pinned
and checked. No native inference, phylogeny or accuracy scoring is performed.

All five canonical traces contain **8,440 constraints**. Their directed
source/target gene sets match the historical trace exactly, including the
constraint sequence order. There are no missing, extra or reordered semantic
constraints. None of the raw JSON trace files is byte-identical to the
historical trace, so raw and semantic equivalence must remain distinct.

The frozen `apply_satellite_membership_constraints` implementation reads only
`source_genes` and `target_genes` from each constraint, converting them to sets.
It tests these sets against high-confidence ortholog pairs from the supplied
family outcomes. It does not consume recorded HMM support, margins, iteration
numbers or cluster labels. Those unused metadata differences therefore cannot
change this consumer's decisions when family outcomes are held fixed.

Constraint order can affect output group ordering in general. A targeted
fixture permutes three constraints with one supported and two detached sources:
membership sets, ortholog pairs, confidence rows and audit counts stay equal,
while ordered group tuples take two forms. The real canonical traces preserve
the historical semantic sequence, so that ordering caveat does not introduce
a difference at this boundary for identical family outcomes.

Thirty-two focused tests pass across the new semantic comparator, existing
phylogeny pipeline tests and residual trace tests. The comparator preserves
direction and multiplicity, rejects invalid/overlapping gene lists, and does
not confuse irrelevant metadata or within-set order with changed constraints.
The test exercises the actual frozen membership consumer, not a replacement
implementation. The source hash matches the installed frozen package.

This closes the candidate-file/constraint-input equivalence checks for the
canonical experiment. It does **not** prove freshly inferred alignments, gene
trees, species trees or pre-filter family outcomes are identical, nor does it
establish final F1 reproduction. No historical score or production default is
changed. Next prepare a bounded fresh phylogeny-stage validation from the
canonical full-fresh candidate output in the frozen installed runtime, without
reusing historical tree checkpoints or rerunning all-to-all search. Pin its
configuration before execution and retain any remaining difference.

```bash
python -m benchmark_tools.audit_canonical_membership_constraints --repo . \
  --output /tmp/canonical-membership-constraints.json
```
