# Bound Post-Measurement Process Review

The single-attempt executor now invokes the process-stream audit after the
measurement wrapper returns and its preflight child has been cleaned up.
It retains the resulting report's path, size and SHA-256 in `result.json`.
A negative sampled-process verdict produces `executor_failed`, preserving both
the wrapper result and negative review. It never retries, submits the next
identity or admits timing evidence. Missing or invalid evidence also fails.

Before creating an attempt, the executor requires explicit finite nonnegative
`maximum_foreign_average_cores` and positive `maximum_sample_period_s` values
in the already bound environmental policy. No real values have been selected;
the test values are synthetic. This adds a prospective sampling-period field
to the policy requirements and must be included in the future reviewed freeze.

The audit binds the following evidence before review and rechecks it afterward:

- Reviewed environmental policy and its v2 process policy.
- Passed preflight response, matching host-policy plan, job, run index and
  exact environmental-policy reference, plus its supporting evidence.
- Collector `ready.json`, native `done.json` and raw `host_processes.jsonl`.

The job scope is derived from the ready worker's actual recorded Slurm
membership; native monotonic times come from `done.json`, not a separately
supplied duration. The typed stream must match the preflight boot, retain the
same observer, bracket native execution and satisfy every adjacent comparison.
Stored interval verdicts are not trusted. Source identity is included, and the
executor's complete recipe checks still run afterward. Reports cannot overwrite
an existing review. Hash binding detects byte changes during review; it is not
an independent attestation that a recorded observation was truthful.

## Validation

217 focused tests pass. New controls cover file-backed positive and negative
reviews, altered policy/preflight/process-policy/supporting bytes, wrong job
or index, raw-log mutation during evaluation, overwrite refusal, and executor
preservation of the review when the sampled process policy fails. The executor
tests substitute synthetic measurement and audit results to test ordering and
failure handling; separate file-backed audit tests exercise the real reviewer.
This is not a native Slurm handoff or full-scale observer test.

The earlier standalone stream-check limitation is now resolved for the normal
measurement-return path. If preparation or collection raises, the existing
failure path remains; partial logs are retained but are not silently promoted
to a complete process review. The review's CPU/process scope does not establish
whole-run executable/configuration stability, quiet pressure/device behavior or
continuous isolation. Those checks, a newly frozen helper recipe, native and
full-scale validation, and a quiet window remain required. No controlled timing,
new scientific result, DGX operation or unrelated process change occurred.
