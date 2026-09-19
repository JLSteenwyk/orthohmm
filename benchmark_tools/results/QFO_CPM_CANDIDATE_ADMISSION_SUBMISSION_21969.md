# CPM Candidate Admission Submission

Array 21969 independently validates corrected candidate preparation array
21967. Both tasks are confirmed PENDING with 2 CPUs/64 GiB and dependencies
afterany:21967 AND aftercorr:21967. The frozen batch requests bizon, two
serial tasks, 4 hours each, no requeue. It does not depend on cancelled
unstarted array 21964.

Executor: `benchmarks/work/publication_qfo_cpm_candidates_admission_v1`,
commit `348454c4d59704725479029e34a11b4fc6143b7f` (pushed).

```bash
sbatch --parsable \
  benchmarks/work/publication_qfo_cpm_candidates_admission_v1/benchmark_tools/results/qfo_cpm_candidates_admission_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_candidates_admission_v1 \
  348454c4d59704725479029e34a11b4fc6143b7f
```

SHA-256 identities:

- `admit_qfo_cpm_candidates.py`:
  `d32d932a8c11f203217db6412fb708c25927af1ffae6a3bb0e122ed8a7e48868`
- `qfo_cpm_candidates_admission_batch_20260919.sh`:
  `d06b85f8f30e93ad480d774510974291dbeef48bd2e647b20b8ddbd419052f2d`

Expected reports: `benchmarks/work/qfo_cpm_candidates_admission_21969_0.json`
and `_1.json`. Logs: `benchmarks/work/qfo_cpm_candidates_admit_21969_0.log`
and `_1.log`. Index 0 is cpm_low, index 1 cpm_high.

Validation requires completed preparation, pinned corrected executor/source,
exact input/helper/context inventory, unchanged expansion parameters and a
single recorded engine invocation. It checks replay authorization and its
fresh reproduction, numeric checkpoint, runtime and complete unique gene
universe; reruns candidate content auditing and independent merge-trace
reconstruction; and rechecks retained file hashes. The independent numeric
decoder's source is recorded. The candidate timing artifact is retained but
not admitted as controlled comparative efficiency.

Tests: 146 focused cases passed in 6.35 seconds, followed by 55 validator
cases in 3.15 seconds after the final decoder-source record addition. Batch
syntax and staged whitespace checks passed. No real candidate result is
yet admitted. This verifies recorded candidate consistency, not independently
rescored search support or biological accuracy. Inferred phylogeny, native
pair validation and official assessment remain required. Existing 11
Dependabot alerts remain unresolved.
