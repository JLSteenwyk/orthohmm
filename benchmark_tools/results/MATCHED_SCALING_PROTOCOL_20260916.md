# Matched Scaling Protocol

Prospective first panel, frozen before new timing outcomes. This addresses the
practical-efficiency requirement, not accuracy tuning or independent validation.
Existing shared-node timings remain descriptive and are not substituted here.

## Inputs and Runs

Use the twelve checksum-pinned complete OrthoBench proteomes from
`orthobench_factorial_prepared_20260916.json` (SHA256
5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382).
Order basenames by SHA256 of `orthohmm-publication-scaling-20260916-v1`, a newline,
and the basename. Use nested prefixes of 4, 8 and 12 proteomes, without removing
proteins or rewriting FASTA content. Record protein and sequence-character counts,
membership and file hashes. Ordering uses neither accuracy nor runtime outcomes.

Run frozen OrthoHMM high-sensitivity, frozen satellite_v2 with inferred phylogeny,
and OrthoFinder 3.1.5 full inference three times each at every size: 27 planned
runs. Rotate the three-method order by repeat plus size index. A run requires a
fresh output directory; do not share inferred search results between methods or
repeats. Record OS page-cache state as uncontrolled; do not flush global caches
on this shared machine. Method settings and native dependency hashes must match
the admitted baseline or be explicitly recorded as a different method version.

OrthoFinder's pre-phylogenetic checkpoint is not a separate complete run. A
checkpoint duration may be reported only if its timestamp and included stages
are observed directly; never subtract unrelated historical timings. Other tools'
historical resource records remain visible with their comparability limitations.

## Resource and Execution Gates

Use identical 32-CPU, 128-GiB limits and a 24-hour per-run wall limit on the same
machine, with no overlapping panel runs. Check requested and actual affinity,
thread settings, CPU model and native runtime. Schedule an exclusive allocation
where available, but also inspect non-Slurm workloads: scheduler exclusivity alone
does not establish a quiet host. Do not terminate unrelated user processes.

Before execution, freeze exact CLI commands, stage boundaries and a tested
resource collector. Record wall time, total CPU time and simultaneous process-tree
memory, with sampling interval and missed-short-process limitations. Retain GNU
time and scheduler records as additional evidence; neither an individual-process
RSS peak nor tracemalloc is the simultaneous total of all native child processes.
If cgroup accounting is available, document its scope and distinguish page cache
from process RSS. Monitor unrelated host activity throughout each run. Classify
contaminated runs explicitly and retain their measurements; replacements, if
necessary, need a documented rule before inspecting comparative speed.

Separate input preparation, tool inference, output conversion and scoring.
End-to-end inference starts before tool database construction and ends only after
native final output completion. Cached/reused stages must not be mislabeled as
end-to-end. Verify complete valid output coverage separately from process exit.
Record failures and timeouts with consumed resources; do not treat them as fast
successful runs or silently omit them from the planned inventory.

## Reporting and Limits

Show every run, median and range for the three repeats, input sizes, failures and
workload flags. Report wall time and memory separately; do not select each tool's
fastest repeat or fit a universal complexity law to three dataset sizes. Paired
comparisons require comparable completed runs, with exclusions made explicit.

This is one nested series: taxon identity and input size co-vary, and full
proteomes preserve their actual duplication and sequence-length distributions.
Results cannot alone establish general scalability or extrapolation to QfO-sized
inputs. Broader dataset coverage and existing resource limitations remain part of
the publication completion audit. Prepared inputs or passing collector tests are
not evidence that the scaling experiments have completed.
