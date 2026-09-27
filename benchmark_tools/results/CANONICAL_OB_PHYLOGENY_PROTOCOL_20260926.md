# Fresh Canonical OrthoBench Phylogeny

## Question and Frozen Inputs

Does a fresh phylogeny-stage run from the validated canonical full-fresh
candidates recover the historical root-group result? This is a reproducibility
diagnostic on development-exposed data, not independent confirmation or tuning.
No historical result will be replaced by this run.

The plan at `benchmarks/work/canonical_ob_phylogeny_20260926/plan.json`
has SHA256 `2ded17c15cf3510e02f1683973edf211dc8b7f7f86037dcf83e121a948dfafdc`.
Its candidate input is the canonical experiment's `fresh_full_self_control`:
54,445 families covering all 251,378 genes in the 12 OrthoBench proteomes.
Candidate bytes match the historical partition. All 8,440 directed membership
constraints have the historical semantic content and sequence order.

Use the pinned installed scientific sources with the private Leiden 0.11
environment, rebuilt MAFFT 7.525 and retained AVX2 FastTree binary. The plan
records the runtime, scientific sources, DendroPy sources, native tools and
input checksums. No package installation or production default change occurs.

Infer the species tree with minimum-variance rooting; retain species-overlap
root groups, positive-paralogy pairs and the high-confidence-pair membership
policy. Recompute alignments and trees; reuse no historical phylogeny checkpoint.
The expensive all-to-all search and validated candidate construction are not
repeated. Labels and scoring are not passed to the native inference stage.

## Preflight and Execution

The retained 16-gene, four-species fixture completed with three candidate/root
groups, two species-tree alignments, one reconciled family, 36 ortholog pairs
and one supported nonempty membership constraint. All checkpoint-hit counters
are zero, including species-tree reuse. Independent structure, sequence,
event/pair and hierarchy readers pass; this does not establish biological truth.
The preparation receipt is
[machine readable](canonical_ob_phylogeny_preparation_20260926.json).
An initial preparation-only duplicate-PATH TypeError was corrected before the
successful fixture plan; no native attempt was retried. Ninety-eight focused
launcher, replay, pipeline, constraint and independent-reader tests pass.

Commit and push this protocol before submitting one Slurm attempt with 32 CPUs,
128 GiB RAM, an eight-hour limit and no requeue. Start the private Python with
`-I`, an explicit clean environment and the pinned plan hash. Do not resume or
retry a native failure automatically. Preserve the failure and scheduler record.
Collect stage resource records, but report shared-host timings descriptively;
this is not a substitute for the dedicated scaling experiment. DGX remains
deferred and is not accessed.

## Admission and Comparison

After terminal completion, recheck plan/runtime/input fingerprints and native
output provenance before independent structure, sequence, event/pair and
hierarchy audits. Require exact input coverage and zero checkpoint reuse.
Score final root groups against the same 70 reference families with the frozen
OrthoBench scorer. Compare partitions, family-level scores and overall F1 to
both historical 74.10607351873405% and fresh installed 73.82156906618316%.
Retain any residual differences and trace them; do not claim score recovery
from candidate agreement alone. No outcome-dependent parameter changes or
replacement of the publication baseline is authorized by this protocol.
