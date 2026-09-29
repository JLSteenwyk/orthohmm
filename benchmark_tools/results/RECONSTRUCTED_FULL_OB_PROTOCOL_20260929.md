# Reconstructed Base: Full OrthoBench Protocol

Freeze this protocol before submission. The
[prepared plan](reconstructed_full_ob_preparation_20260929.json) is a single
full-data reproduction attempt, not a new accuracy configuration or timing
repeat. Existing full job 22337 validated the recovery workflow with the
original base interpreter. The reconstructed base has only completed the
16-gene installation fixture. This experiment tests that specific remaining
dataset-scale execution path.

## Fixed Scope

- All 251,378 genes in the same 12 proteomes, 70 RefOGs and 11 uncertain files.
- Same frozen scientific package, canonical-order harness, recovery and reader
  locks, MAFFT/FastTree assets, inputs, seeds and 32-CPU settings as job 22337.
- Controller, venv base and offline pip installer use the reconstructed
  `/tmp/orthohmm-base-reconstruction-20260927/python-runtime/bin/python`.
- Fresh inference and reader venvs, private input copies and native outputs;
  no search, clustering or phylogeny checkpoint reuse and no automatic retry.
- One local Slurm job: 32 CPUs, 128 GiB, 24 hours on `bizon`. Shared-host
  correctness only. Do not stop other work or claim controlled resources.
- Run root: `benchmarks/work/publication_reconstructed_full_ob_20260929`.

The plan SHA256 is
`cea0dd8efa459a01c57ca0450798005c5de3aeef15a6fda19c181acd35918bee`.
It checks 6,515 file records before/after execution, including 6,298 regular
base-runtime file payloads (resolved symlink payloads included, pyc/pyo caches
excluded). Original 215 records remain pinned. This inventory is not a proof
of all OS dependencies, symlink topology or unobserved imports. No existing
runtime or historical plan is modified.

Exactly four command arguments change: controller interpreter, base Python,
installer Python and output directory. The timing runtime has different
scientific dependency versions and is deliberately not substituted here.

## Admission And Comparison

Require successful terminal scheduler state, all eight stages, exact planned
commands/environments, unchanged pinned inputs, fresh venv installation
inventories and wheel payload audits, and all four independent scientific
readers. Recompute the complete 70-family score object and compare the entire
root partition without relying on group labels against admitted job 22337.
Compare native pair/event/hierarchy outputs as well; report differences rather
than selecting a favorable run. No parameter tuning or default change follows.

The existing `admit_integrated_full_ob` command is bound to job 22337 and
must not be reused with a substituted job ID. A separate admission must bind
this new plan and eventual scheduler job; controller success alone is not
admission. Preserve failure logs and stop on failure without retry.

Reference expectations are 59,770 root groups; F1 74.10607351873405%,
precision 81.77045380181866%, recall 67.75533630827974%; 966,439 native pairs,
8,681 reconciled and 45,764 bypassed candidate families. Equality is the
reproduction endpoint, not a reason to conceal a changed result.

## Limits

This is not an independent biological test, complete archive restoration,
all-method validation or cross-host experiment. Existing native assets/readers
remain separately supplied. GNU-time output covers the workflow on a busy
host; maximum process RSS is not simultaneous process-tree peak memory.
No result enters the frozen 27-run controlled timing panel. Its isolation and
accounting gates remain unchanged, as do QfO uncertainty and release-rights
limitations.
