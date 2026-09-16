# Saved-graph CPU-affinity diagnostic

## Question and scope

Does CPU availability explain the different initial Leiden partitions observed
with the same preserved QfO RBNH graph? The prior single-CPU instrumented repeats
reproduced one partition three times. Historical runs requesting 32 CPUs produced
a different partition, but their actual worker affinity was not recorded.
Affinity is a hypothesis, not an established cause. No accuracy is evaluated.

## Frozen design

Use the admitted capture from job21305: identical gene order, integer endpoints,
float64 weights, CPM resolution0.1, seed4, and include-isolates enabled. Preserve
the existing frozen49ab110 worker and native runtime. The observer imports and
records the same modules and native libraries in every arm.

Within one32-CPU/64GiB allocation, run four fresh workers sequentially:

1. `one_cpu_0`: lowest allocated CPU.
2. `all_cpus_0`: lowest32 allocated CPUs.
3. `one_cpu_1`: same single CPU as arm1.
4. `all_cpus_1`: same32 CPUs as arm2.

Reject allocations with fewer than32 distinct available CPUs. Each child sets
its own affinity before importing the frozen worker or native scientific
libraries. Record inherited, requested, and actual affinity; reject mismatches.
Fix PYTHONHASHSEED=0 and OMP/OPENBLAS/MKL thread variables=1 in every arm.
Changing affinity does not request32 BLAS/OpenMP threads. Parent affinity stays
unchanged. Failures remain preserved without automatic retries.

Hash-check graph inputs, source modules, loaded libraries, interpreter and frozen
runtime. Require equal recorded software identity excluding the intentionally
changed affinity. Retain all partitions and compare complete membership against
the capture, diagnostic, first repeat, and preceding same-affinity arm.
Do not choose a result based on benchmark performance or promote new defaults.

## Interpretation limits

Matching partitions within each affinity but different partitions between arms
would implicate CPU availability under this instrumented context; it would not
establish the internal mechanism or reproduce unrecorded historical binaries.
Identical partitions across arms would leave process/import context and other
unmeasured historical differences open. Nonrepeatability within an arm requires
further diagnosis before declaring a reproducible publication baseline.

Four alternating runs are a bounded diagnostic, not proof of general determinism.
The shared node has unrelated workloads, so elapsed times are not controlled
efficiency benchmarks. No graph rebuilding, profiles, inference accuracy, or
parameter optimization is included. Full QfO ablations remain pending.
