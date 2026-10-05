# Remaining Native Request Handoff

`prepare_native_factorial_request.py` replaces repeated per-index preparation
scripts for indices3-12 without changing the frozen executor, plan, policy,
scientific method or any existing request. It never submits, releases,
cancels, retries or resumes a job. Calling it is not a native-start claim.

The index is the length of the explicit reviewed prefix, not an arbitrary
requested cell. Each pinned history entry is checked through the existing
`reviewed_history` validator, including the explicit failed-index-zero
adoption and fresh Slurm terminal outcomes. Require a new held job ID larger
than every previous ID, the exact shared32-native-core/128GiB envelope and
absent next-run/session outputs. Unsafe available RAM refuses preparation;
ordinary background activity remains diagnostic-only.

For a preceding successful OrthoBench identity, also require its separately
bound score. QfO native output validation permits the next inference identity
without inventing a QfO score: pair conversion/native six-endpoint assessment
and independent admission remain separate scientific work. Failed attempts
retain their failure; this helper adds no retry or replacement authorization.

Use exact original review references from the running request history. A
published byte-identical copy is not automatically an identical pathname
binding for its separately pinned score. Do not rewrite old receipt paths.

After the preceding identity is actually terminal, fully reviewed and scored
where required, submit only its successor on hold using the unchanged batch:

```bash
sbatch --hold --parsable --partition=gpu \
  --output=benchmarks/work/native_factorial_launch_20261004/index03_%j.out \
  --error=benchmarks/work/native_factorial_launch_20261004/index03_%j.err \
  benchmark_tools/run_native_factorial_cost.sh \
  --request /ABSOLUTE/FRESH/request_03_receipt_amended.json \
  --request-sha256 scheduler-comment
```

Prepare the request under the retained sanitized Python3.10 environment:

```bash
PYTHON -B benchmark_tools/prepare_native_factorial_request.py \
  --job NEW_HELD_JOB \
  --history ORIGINAL_INDEX0_REVIEW SHA256 \
  --history ORIGINAL_INDEX1_REVIEW SHA256 \
  --history ORIGINAL_INDEX2_REVIEW SHA256 \
  --prior-score ORIGINAL_INDEX2_SCORE SHA256 \
  --output /ABSOLUTE/FRESH/request_03_receipt_amended.json
```

The preparation pins its own source and current Git revision, fixed amendment,
full history, score if required, held scheduler observation and available-memory
snapshot. Its creation-only output retains evidence without altering prior
requests. Tests use synthetic scheduler responses; they do not establish a
real held-job handoff. Actual preparation must pass these checks on the live
host before release.

Bind the resulting digest in the owned held job's Comment, then independently
reread the request, full prefix, held scheduler identity and safe capacity
before releasing that one job. Verify actual controller/native process and
resource limits before reporting inference started. This helper intentionally
does not combine those actions. Never release based on a startup receipt or
an unreviewed apparent output. The native controller still repeats its own
launch/runtime/resource checks; no existing gate is weakened.

These are shared-host observations with unknown, potentially tool-dependent
timing distortion. No quiet window, DGX, isolated-speed claim, background
subtraction or change to unrelated scientific jobs/services is introduced.
