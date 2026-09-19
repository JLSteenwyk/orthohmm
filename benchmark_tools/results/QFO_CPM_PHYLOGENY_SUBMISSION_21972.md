# CPM Inferred-Phylogeny Submission

Array 21972 runs cpm_low (task 0) and cpm_high (task 1) serially. Both tasks
are confirmed PENDING with 32 CPUs/192 GiB and dependencies afterany:21969
AND aftercorr:21969. The frozen batch specifies bizon, 24 hours, no requeue.
It waits for the entire admission array to terminate and corresponding
candidate admission to succeed; failures must not be bypassed.

Executor: `benchmarks/work/publication_qfo_cpm_phylogeny_v1`, commit
`f1a2bb875f5a56a827e5669780b2f1c8718dd1cd` (pushed).

```bash
sbatch --parsable \
  benchmarks/work/publication_qfo_cpm_phylogeny_v1/benchmark_tools/results/qfo_cpm_phylogeny_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_phylogeny_v1 \
  f1a2bb875f5a56a827e5669780b2f1c8718dd1cd
```

SHA-256 identities:

- `run_qfo_cpm_phylogeny.py`:
  `f33e02be25e91728ea83edcf64eb021b990e324a9a9036657ca7fbc2511b63d6`
- Shared baseline verifier `run_qfo_parameter_phylogeny.py`:
  `001afb28ac671d454883f3ff80c9424dc474596d652cfad9b0ce73dee77db779`
- `qfo_cpm_phylogeny_batch_20260919.sh`:
  `ea1249099e721495f28ebf3a90048f95f2e52b8cbbe8b0d1b410894e236e346b`

Outputs: `benchmarks/results/qfo_cpm_phylogeny_v1/cpm_low` and `cpm_high`.
Scheduler logs: `benchmarks/work/qfo_cpm_phylogeny_21972_0.log` and `_1.log`.
Each output retains preflight, a freshly reproduced independent candidate
admission, native execution evidence, and successful or failed postflight.

Both cells use their own admitted CPM seed/candidate/constraint artifacts.
The species tree is inferred independently, not supplied from the baseline.
The unchanged native engine may reuse raw alignment/tree checkpoints only
under its exact membership, sequence and tool-identity rules. Candidate
admission source/commit and transitive file records are checked before
launch; baseline native artifacts, command-source equivalence and runtime
are verified using the existing baseline checks. Source and input identities
are checked again after execution. Failures are retained without retry.

Validation: 109 focused tests passed in 1.03 seconds; batch syntax and staged
whitespace checks passed. This is submission, not native-output validation
or accuracy admission. Native reconciliation/pair validation, conversion
and official scoring remain required. Shared-host cached runtime is not
controlled efficiency evidence. No scientific default or publication claim
changed. Existing 11 Dependabot alerts remain unresolved.
