# Prospective Threadripper Timing Amendment

## Scope

The user selected the existing local Threadripper instead of the DGX. This
amends the host and execution resources of the replacement timing panel;
it does not change scientific settings, inputs, endpoints or publication scope.
Do not access the DGX or make its service permissions a prerequisite.
Historical ARM and shared-host measurements remain separate descriptive data.

Preserve all 27 identities and their frozen order: high-sensitivity OrthoHMM,
satellite_v2 OrthoHMM and full OrthoFinder 3.1.5, each with 4/8/12 complete
proteomes and three repeats. No protein subsampling, search reuse, output
reuse, best-repeat selection or changed scientific defaults is allowed.
The [input parity receipt](threadripper_input_parity_20260928.json) verifies
the 27 identity/order tuples against both retained command plans and all 12
local input-file hashes against the DGX input manifest. The subsets contain
73,266, 165,168 and 251,378 proteins. Native input enumeration must also match
the retained order; filename-set equality alone does not establish this.

## Allocation And Placement

Use node `bizon`, partition `gpu` with no GPU request, one sequential exclusive
node allocation per run, 32 native worker CPUs and 128 GiB memory, a 24-hour
Slurm limit and 23h50m native timeout. Disable automatic requeue and retries.
Do not modify global Slurm settings. Verify actual job allocation and enforced
cgroup limits before any inference, not merely requested submission flags.

The host has a Ryzen Threadripper PRO 7995WX, 96 physical cores, 192 logical
CPUs, one NUMA node and 1,081,211,752,448 bytes of RAM. The full `scontrol`
record reports RealMemory=1030000 MB, so 128 GiB is feasible. The narrower
`sinfo -N -l` display truncates that value to 103000; do not infer a smaller
memory capacity from this display. Current Slurm topology reports 192 sockets
with one core each, unlike physical topology. Do not use those scheduler
socket labels as physical-core evidence.

Bind native inference and descendants to the same 32 logical CPU IDs drawn
from distinct physical cores, validated against Linux topology and the
allocation's effective cpuset. Record the chosen IDs before the first native
run; keep them unchanged across the panel. If scheduler binding cannot provide
that placement or cgroup limits cannot be verified, do not launch. Retain OS,
kernel, CPU-governor, tool/runtime and affinity identities for comparisons.
Do not equate 32 allowed CPUs with proof of 32 continuously utilized cores.

## Quiet-Window Policy

No concurrent unrelated scientific analysis is permitted for controlled
timing, even if it is outside Slurm or uses different CPU IDs. Check the host
process inventory as well as scheduler state before each run. Freeze the
ordinary background-service inventory before the first run. Do not silently
stop, suspend or alter unrelated processes or services. Request coordination
when those workloads prevent a quiet window.

Bracket each native run with whole-host observations and continue observations
through execution. Retain raw process identity/cgroup CPU deltas, observation
gaps and errors, pressure and available memory/I/O evidence. Monitoring is
not proof that every short-lived process or source of contention was seen.
Known overlapping scientific work, configuration/placement drift or missing
required monitoring evidence makes a run ineligible for controlled comparison
and stops subsequent submissions pending review. Do not hide its measured
resources. Ordinary native CPU or memory stalls are not contamination criteria.
Do not subtract estimated foreign workload or collector overhead from times.

The existing host-CPU detector remains a diagnostic, not an automatic quiet-host
certificate. Its 0.25-core threshold is not a new threshold for deleting slow
runs. Resolve background-workload eligibility independently of method outcomes.

## Measurement And Failure Handling

### JIT Cache Policy (2026-09-28)

Before production timing, use a fresh per-run `NUMBA_CACHE_DIR` beneath the
persistent run directory for every method. Create it empty outside inference;
include any native JIT compilation and cache writes inside inference time.
Use different fresh directories for before/after verification probes so they
cannot populate the native cache. Retain resulting cache-file inventories;
do not remove or reuse historical shared caches. This is cold JIT-cache
execution on prepared, hash-warmed FASTA inputs, not cold filesystem caching.
Full runtime checks still require isolated Python bytecode lookup separately.
Restrict `NUMBA_CACHE_LOCATOR_CLASSES` to `UserProvidedCacheLocator`, supported
by the installed runtime, so an unusable explicit directory cannot fall back
to a historical source-tree or user-wide cache.

This closes a prospective gap found before any production timing: the frozen
source has `@njit(cache=True)` and the shared source directory already contains
`.nbi`/`.nbc` files. Python's bytecode controls do not define Numba cache state.
See [Numba's environment reference](https://numba.readthedocs.io/en/stable/reference/envvars.html#numba-cache-dir).
Earlier fixtures remain valid integration controls, not cold-JIT timing evidence.
The native collector, CPU settings, scientific algorithms and run order are
unchanged. Validate this policy on the bounded fixture before production.

### Input Storage Clarification (2026-09-28)

Before any scientific timing result, adopt one common preparation policy:
fresh per-run FASTA copies on local `/dev/shm` tmpfs for all three methods,
with results, search databases, native intermediate files and collector
artifacts on the existing project ext4 volume. Create and hash the input
copies outside the inference timer, in the scheduled job's preparation
context, not in an unrelated long-lived process. Do not use the exploratory
tmpfs copies directly as production inputs or share a mutable input directory
between runs. Verify actual native enumeration before and after every run;
creation order alone is not sufficient. Preparation time and copied bytes
are separately reported. This is a prepared-input measurement, not cold-disk
end-to-end performance; do not evict system caches.

Use explicit OrthoFinder `-o <persistent-run-directory>/orthofinder_full`
alongside the fresh tmpfs `-f` input. This changes output placement only;
retain all frozen scientific flags, including full phylogenetic inference.
Do not use `-op` in production: it is solely the preparation-only diagnostic
used to verify initial path behavior. That control passed on 28 September
and left its two synthetic input files unchanged. Its
[receipt](orthofinder_storage_probe_20260928.json) pins native outputs and
commands. Full-run temporary-file placement remains a validation gate.

The job-level 128-GiB cap must include preparation and native work. Record
job and native-step memory scopes separately. Pre-existing tmpfs pages may
remain charged to their creator rather than the reader, so native-step peak
alone is not the complete memory footprint. Report input bytes and job peak
alongside native-step peak, without claiming their peaks are additive or
subtracting prepared input memory. This policy requires an accounting probe
and executor integration before production; the earlier sleep controls do
not verify tmpfs charge ownership. Inputs are temporary, not archival data.

Record input preparation, copying/hashing, inference, output validation,
conversion and scoring as separate stages. Native timing spans CLI launch
through exit, including database construction, descendants and output writing.
Use the same tested cgroup/process-tree accounting for all methods, reporting
wall time, actual CPU time, peak memory and their accounting semantics.
Do not substitute GNU-time maximum single-process RSS for aggregate tree memory.
Collector validation on this x86 host is required before production runs;
previous ARM collector checks alone do not establish it.

Preserve terminal status and consumed resources for every failure or timeout.
Observe a live job until authoritative terminal evidence is available; an
observation timeout does not authorize a duplicate run. Native failures may
advance to the next identity only after terminal/output checks and continued
environmental eligibility. Infrastructure or policy failures stop submissions.
No replacement of an unfavorable observation is authorized by this amendment.
Report incomplete repeats as incomplete, never as complete three-run summaries.

## Current State And Next Execution Gates

The three-method native installation fixture completed in job 22350 after
two retained preflight failures. The local
[baseline amendment](threadripper_native_baseline_20260928.json) corrects only
cwd-dependent OrthoFinder package metadata: historical inventory included
repository-local `orthohmm.egg-info`; actual frozen execution cwd does not.
No installed package or scientific source was changed. Launchers must explicitly
unset inherited `LD_LIBRARY_PATH`, `LD_PRELOAD` and `LD_AUDIT`; runtime checks
still reject loader overrides rather than silently accepting them. The
[fixture receipt](threadripper_native_fixture_22350.json) records successful
native output checks and independent resource replay, but no timing admission.
Its whole-host samples measured 72.71-73.67 competing CPU-core equivalents.
Full transitive runtime checks and production execution gates remain open.

A read-only observation on 28 September found an empty Slurm queue but about
69.6107 CPU-core equivalents of persistent outside work. IQ-TREE, BAli-Phy,
HyPhy and Python analyses are active under user services and interactive
scopes. No process was signalled and no native benchmark was launched.
The raw preflight stays local at
`benchmarks/work/threadripper_preflight_20260928.json` (1466021 bytes), SHA256
`b217285f111b4671f7ad44efab440e87191d7a326d58696d46e13bd110c3003a`.
Its CPU observation excludes only the observer, uses per-process intervals,
and had zero sampling errors; short-lived processes may be missed. It is not
a whole-run measurement or a prediction of future load.

Next: validate the retained x86 runtime and deterministic enumeration; deploy
the local executor and verify allocation/binding/collector behavior; secure a
quiet window; then run the frozen complete panel. This amendment freezes
prospective resource and failure policies, not execution authorization under
the currently busy host conditions. Publication readiness remains unproven.
