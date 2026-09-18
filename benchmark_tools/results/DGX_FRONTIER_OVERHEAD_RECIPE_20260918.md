# Verified Native Overhead Recipe

Launcher committed/pushed as75b7723. The new remote recipe is isolated at
`/home/jlsteenwyk/projects/orthohmm-publication/frontier_overhead_recipe_v1`.
Its39 files match local sources byte-for-byte; no external symlinks are
present. Runtime/source checks remain required before and after every native
command, outside its measured wall interval.

Recipe manifest SHA-256:
`64a05b5201e78a9d8d46f302879a49bd5cc470bf7813eb5d16e4c88d79c51edc`.
Separate exact-scope authorization SHA-256:
`77875e0273883454c25fcb697916aedc10f5d0d3081baa77ead997d6f4aa0b47`.
Plan SHA-256:
`58e482e6f123bc82de5e98df3c93f4bf0ca3a47d2302ffcaa47290c7d1e20685`.

Authorization covers indices0..17 of this incremental-observer-overhead
experiment only, with scientific_execution_authorized=false. It binds the
plan and transferred recipe exactly. The original plan's execution flag
remains false; the explicit separate record enables only this experiment.
The launcher additionally checks that itself, both collectors and the plan
are present with matching hashes in the pinned recipe. It retains per-task
authorization/mode provenance separately from native measurements.

A DGX preflight imported the transferred launcher and successfully selected
all18 tasks using these exact pins, with native_execution_started=false.
The local scheduler reported no DGX jobs before submission preparation;
the filesystem had3.2TiB available. These observations are not a guarantee
of absence of unscheduled activity. No unrelated services/jobs were stopped.

58 focused tests pass, covering all18 selections, strict authorization
types/scope, recipe-bound collectors, both observation modes, failed native
status preservation, paired-command derivation and boundary lifecycle.
Shell syntax validation passes. No overhead measurement is claimed yet.
