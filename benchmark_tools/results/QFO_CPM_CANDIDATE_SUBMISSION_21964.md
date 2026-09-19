# CPM Candidate Preparation Submission

## Superseded Before Execution

Array 21964 was administratively cancelled before either task started.
Accounting confirms both tasks CANCELLED by 1000, Start=None, zero allocated
CPUs and 00:00:00 elapsed. Code inspection found the numeric auditor imported
the executor's scientific package before frozen launcher selection. Its
source guard would reject execution; the mocked tests had missed the import
order. The corrected executor delays that auditor import and adds a real
fresh-process regression check. No result-based retry, native output or
score selection occurred. Original submission details are retained below.

## Original Submission

Array 21964 has two serial tasks: index 0 cpm_low, index 1 cpm_high.
Scheduler inspection confirms both PENDING, 2 CPUs/64 GiB, with
afterany:21962 AND aftercorr:21962 dependencies. Frozen batch specifies
bizon, 4 hours, no requeue. The whole upstream admission array must be
terminal and the corresponding task successful. Do not bypass failed
dependencies or rerun selectively based on outcomes.

Frozen executor: `benchmarks/work/publication_qfo_cpm_candidates_v1`,
commit `6e209994961ac71483603a28a73bb0d56c91d013` (pushed to origin/main).

```bash
sbatch --parsable \
  benchmarks/work/publication_qfo_cpm_candidates_v1/benchmark_tools/results/qfo_cpm_candidates_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_candidates_v1 \
  6e209994961ac71483603a28a73bb0d56c91d013
```

SHA-256 identities:

- `prepare_qfo_cpm_candidates.py`:
  `a87a291e3987983f3e887c01516a0f78cb5940fe0ea377e8f7901f34cf35dab3`
- `qfo_cpm_candidates_batch_20260919.sh`:
  `b296891b971b0754892961705278087443b64119d3cf13c24bb6888fb1ef6dbd`

Each task requires its admitted CPM-specific strict_profiles_refined seed,
not the baseline seed. The existing frozen satellite_v2 expansion uses
unchanged candidate parameters; the CPM change is the only perturbation.
The independent replay admission is freshly reproduced before construction.
The existing candidate auditor checks seed membership, complete partition
coverage, sidecars and merge-trace reconstruction. Runtime/checkpoint/input
identities are checked before and after. Failures preserve partial outputs.

Outputs: `benchmarks/results/qfo_cpm_candidates_v1/cpm_low/manifest.json`
and `cpm_high/manifest.json`, with scientific artifacts under each
`candidate/orthohmm_working_res/` directory. Logs and GNU time records are
`benchmarks/work/qfo_cpm_candidates_21964_{0,1}.log` and `.time.txt`.
The time records include verification overhead on a shared host; they are
not controlled end-to-end inference measurements.

Validation: 128 focused tests passed in 5.37 seconds; batch syntax and
staged whitespace checks passed. Tests include variant-specific seed use,
unchanged parameters, content auditing, source/context/scheduler gates,
fresh-admission mismatch and retained failure reports. No real candidate
result is admitted yet. Independent candidate admission, independently
inferred phylogeny, native pair conversion and benchmark scoring remain
required. No default or publication claim changed. Existing 11 Dependabot
alerts remain unresolved.
