# Relocated OrthoBench Factorial Reproduction

## Executed Result

Exported nine files from commit
`a8a9f57734902715df7e9220cba2b3e0106fd5f3` into a separate directory and
executed the exported analysis with isolated Python (`-I -B`). The run
completed with exit code zero and empty stdout/stderr. All scientific fields
exactly match the committed reference result, including all eight cells,
paired effects, confidence intervals and interaction statistics.

The analysis uses retained sufficient statistics for 70 RefOGs, 20,000
shared paired draws and seed 20260918. This remains development-exposed
evidence. No parameters, endpoints, scores or interpretation were changed.

Machine-readable evidence: `orthobench_factorial_relocated_20260918.json`.
It records the exported source/data identities, interpreter and packages,
command, output hashes and limitations. All nine exports and both outputs
were rehashed after completion; a separate comparison against the retained
worktree result also passed. Nineteen focused reproduction/bootstrap tests
pass.

## Reproduction Command

From the repository root, use an analysis virtual environment with the
versions in `benchmark_tools/swiss_analysis_requirements.txt` and the
corresponding hash lock. The executed environment used Python 3.10.13 and
NumPy 2.2.6. The runner records installed versions; it does not install or
enforce the dependency lock. Choose new output/report paths on every run.

```bash
python benchmark_tools/reproduce_orthobench_factorial.py \
  --revision a8a9f57734902715df7e9220cba2b3e0106fd5f3 \
  --python benchmarks/work/swiss_analysis_env_20260917/bin/python \
  --output /tmp/orthohmm-ob-factorial-reproduction-new \
  --report /tmp/orthohmm-ob-factorial-reproduction-new.json
```

The original execution output is
`/tmp/orthohmm-ob-factorial-reproduction-20260918/reproduced/`.
It contains `statistics.json` and `statistics.md`. These temporary files
can be regenerated from the committed inputs; their presence is not a
long-term archival guarantee.

## Scope

This demonstrates relocation of the statistical workflow without loading
the historical absolute input paths stored as provenance. It does not
rerun native inference, partition conversion, reference construction or
official scoring. It does not establish cross-platform reproducibility,
an OS-hermetic environment, reference-data redistribution rights or full
publication readiness. Those requirements remain separate.
