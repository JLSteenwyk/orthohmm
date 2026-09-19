# GPU-Available Long-Target Routing Fix

## Defect And Fix

Source inspection during preparation of search-rejection tracing found an
uncovered path in `search_species_pair_indexed`: CUDA could be available
while every candidate target exceeded its1998-residue limit. In that case
`n_gpu == 0`, but the CPU-only branch also required `not use_gpu`. Neither
backend ran, and the uninitialized score array reached E-value estimation.
This is a correctness defect, not a proposed accuracy optimization.

The current engine now runs CPU scoring whenever a nonempty batch has no
GPU-scored pairs. A failed GPU call also reaches that branch exactly once;
mixed successful GPU/CPU batches keep their original routing. Kernel code,
score formulas, candidate filtering and default parameters are unchanged.

## Regression Evidence

Before the fix, the new12-case routing test had2failures and10passes. Both
failures had CUDA availability true and only1999-residue targets. The test
also checks backend call inventories, so it cannot pass accidentally just
because uninitialized memory happens to match the expected score.

After the fix, the same routing tests plus existing engine tests passed28/28.
Expanded coverage exercises multipair C, scalar C fallback and Numba paths,
CUDA unavailable/available/failure,1998/1999boundary targets, mixed batches
and target-order reversal. A further real CPU-scoring fixture compares
scores and E-values for1999/2001-residue targets with CUDA availability
forced off versus on. Current C scoring availability was confirmed true.
The final focused run passed53tests in1.17seconds:

```sh
python -m pytest -q tests/unit/test_engine_score_routing.py tests/unit/test_engine.py \
  --junitxml=benchmarks/work/search_routing_native_cpu_20260918.xml
```

Retained JUnit evidence in `benchmarks/work/`, SHA-256:

| Report | SHA-256 |
| --- | --- |
| `search_routing_before_20260918.xml` | `be7c16ecc90d3e52dc45262bd61e170fa227282c81fb861835bc494043a54134` |
| `search_routing_after_20260918.xml` | `e46e26e9c8faea007b50d5039b7fac8c1385fd32febe263693a59c4c68a02cc2` |
| `search_routing_all_backends_20260918.xml` | `92ecfb8c204c1ce4823d80db9093eb0836480e27c25b9825daae21ff7e2fdc0c` |
| `search_routing_native_cpu_20260918.xml` | `a65d0202ab2660db2cb202b66ad10b233d967b93c83c6f725fc93e623b62e4e0` |

GPU routing tests mock GPU execution; this is not numerical GPU/CPU kernel
parity evidence or a full-suite rerun after the fix.

## Scientific Impact Boundary

The frozen publication core also contains the defective branch, but its
[retained native runtime](publication_native_runtime_20260916.json) declares
`cuda_enabled: false` and exactly three CPU libraries. A fresh runtime
verification rehashed all recorded sources/binaries and checked that exact
library inventory; no CUDA library is present. The frozen loader returns
CUDA unavailable when its package-local CUDA library is absent. The corrected
QfO command plan binds this CPU-only runtime, and its admission checks that
runtime. This evidence supports exclusion of this trigger for the admitted
CPU-only configuration; it is not a historical process-level execution trace.

Do not extend that conclusion automatically to older GPU-enabled benchmarks,
installed wheels, other runtime directories or DGX jobs. They require their
own runtime and candidate-batch evidence. No earlier accuracy score is changed
or declared invalid here, and no broad historical non-impact claim is made.
Frozen executors, active jobs, native libraries and reference data were not
modified. This fix is for the current source tree, not a silent replacement
of the evaluated scientific core. Any new scientific configuration still
needs explicit freezing and applicable independent confirmation.

The planned prefilter-versus-scoring rejection trace remains unfinished;
this defect is not assumed to explain the retained OrthoBench misses.
