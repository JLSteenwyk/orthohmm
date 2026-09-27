# Fresh Retained-Order QfO Phylogeny Preparation

Implements the first native arm of the previously frozen
[downstream comparison](QFO_ORDER_DEPENDENCY_TRACE_20260927.md).
The [preparation receipt](qfo_fresh_phylogeny_preparation_20260927.json)
binds the full plan and small-fixture scientific checks.

Full plan: `benchmarks/work/qfo_fresh_phylogeny_20260927/plan.json`,
1,035,097 bytes, SHA-256
`57247bf76d26c1ff54cf13aa4c9d0a9622f0b8a2902c45ef222842513888d030`.
It pins admitted retained-order candidates/constraints, all 78 proteomes,
installed recovery package payloads, historical comparison products and
validated MAFFT/FastTree executables/helpers. Same frozen scientific rules;
32 local CPUs, 128 GiB, six hours, one attempt and no requeue/retry.
The isolated worker imports the installed package before repository helpers.
It records native arguments/environment and all output-file identities.

Validation before submission: 42 focused tests pass. The 16-gene four-species
fixture completed fresh inference: three candidate/root groups, one reconciled
family, two bypassed families, 36 native pairs, four duplications, three
speciations, two inferred-tree markers and zero checkpoint reuse. Structural,
sequence, event and hierarchy readers pass. Empty fixture constraints do not
validate the full QfO merge trace or establish full-data equivalence.

The full arm is not submitted at this preparation milestone. Commit/push
before submission; record the real scheduler ID afterward. Require successful
scheduler termination, unchanged artifacts, full scientific readback and
historical output comparison before admitting any canonical reuse. The
canonical launcher/admission remains to be implemented separately; neither
fresh-arm completion nor accuracy equivalence is claimed here.

This reuses expensive admitted search/candidate results, not historical trees.
Timing is descriptive on the shared host. The DGX remains deferred and the
publication goal is incomplete.
