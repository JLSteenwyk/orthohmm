# Full Recovery OrthoBench Protocol

One fresh full OrthoBench run will test the explicit canonical-order pipeline
entrypoint in the newly installed recovery environment. The
[prespecified protocol](full_recovery_orthobench_protocol_20260926.json) binds
the prepared plan, original reference files and historical baseline. This is
a development-exposed reproducibility experiment, not independent validation
or an accuracy improvement search. No settings will be tuned from its result.

Inputs are the verified 12 proteomes containing 251,378 unique proteins. Search
starts from FASTA; no search, grouping or tree checkpoints are reused. The
installed frozen scientific package uses Leiden 0.11, built-in HMM search,
high-sensitivity refinement, CPM resolution 0.1, satellite-v2 expansion,
inferred/minimum-variance species rooting, species-overlap gene rooting and
positive-paralogy pairs. Only candidate input ordering is explicitly adapted
by the separately recorded canonical policy. Production defaults do not change.

The launcher requests 32 CPUs, 128 GiB RAM and 24 hours on the shared host,
with no requeue. It allows one attempt, rejects existing native output, records
failures, and does not automatically retry or resume. It checks the plan digest,
executor, all existing top-level benchmark Python sources, installed wheel
payload, source/FASTA/tool records and dependency audit. The native command
runs the isolated recovery interpreter. The complete installation audit is
repeated before and after inference; GNU time and native logs are retained.
These resources are descriptive, not controlled comparative timing. DGX is
not accessed.

## Readback Gates

Require terminal scheduler success, executor completion and independent native
artifact binding. Run structural, sequence, event/pair and hierarchy readers.
Compare the root partition with the historical baseline, ignoring group IDs;
recompute the complete frozen score object and all 70 RefOG records using the
81 pinned reference/low-certainty files. Preserve all discrepancies and
explain them without changing the endpoint or replacing earlier scores.

The current plan is at
`benchmarks/work/publication_full_recovery_orthobench_20260926/plan.json`, SHA256
`57878c0712b3216376ea1d431b1d8e9e8e3ae8eb58e3bd3069e15cf7dddb7005`.
The plan and its referenced sources must not change during the run. New
readback files may be added separately. Twenty-seven focused launcher,
ordering-wrapper and package-audit tests pass before submission. The previous
full-pipeline fixture passed 86 focused tests and all four scientific readers,
but made no satellite merges; that result alone does not validate this run.

No new scores are admitted by preparing or launching this experiment. Broader
generalization, uncertainty, timing and publication-package requirements remain
open even if exact baseline recovery succeeds.
