# Environmental Handoff Worker

`benchmark_tools/threadripper_environment_worker.py` implements the separate
preflight responder required by the local executor. It never submits a job,
writes the native `go.json`, retries a failure or authorizes the next identity.
This is implemented and component-tested, **not production-validated**.
No real policy, readiness approval, execution permit or native run was created.

## Inputs And Operation

Run the worker within the same local Slurm allocation as the executor, before
the executor reaches its parked-worker handoff. It takes the execution request
and environment policy paths plus their SHA-256 digests. The existing selector
checks the frozen command/lookup plan, complete helper-source recipe, reviewed
history position and readiness. The policy must also be bound directly by
`readiness_review.environment_policy`; a matching boolean is insufficient.

The worker waits up to 4,000 seconds for the measurement directory's release
request. The native guard's existing 20-second response deadline starts when
that request is written, not when this waiting worker starts. Wrong identities,
stale requests, preexisting responses and source drift cannot authorize release.
Waiting or validation failure does not launch a replacement attempt.

The policy uses schema `threadripper_environment_policy_v1` and binds:

- The local host and frozen plan, `decision: reviewed`, an explicit review
  reference and supporting evidence records.
- A v2 process-policy file with same-boot, reviewed process identities. Each
  ordinary user-process entry additionally has an `image` file record for its
  approved executable. Kernel entries retain the narrow typed-name rule.
- Nonempty `configuration_files` for the independently reviewed ordinary
  service/configuration inventory. File hashes do not prove that the inventory
  is complete or that its classification is correct.
- `maximum_foreign_average_cores` and `maximum_pressure_percent` for CPU, I/O
  and memory. These are prospective policy values, not defaults selected by
  this implementation. The 0.25-core diagnostic threshold is not reused as an
  admission rule. Real values still need to be fixed and reviewed before timing.

At handoff, the worker verifies the collector's initial flushed process sample,
live observer and parked native identities, actual job/step memberships, and
absence of a native release or completion file. It collects two typed process
snapshots separated by three seconds and retains host pressure/counter reads.
It applies the reviewed process policy, the explicit background CPU bound and
the pressure bounds using the existing PSI parser and interval calculation.
Any unresolved process, missing type, sampling error or exceeded bound fails.

For each approved ordinary process, it hashes the executable reached through
`/proc/PID/exe`, checks it against the approved image bytes, and rechecks process
identity and executable inode metadata. Shared identical loaded images reuse
the digest within this check; replacement, metadata changes or PID drift fail.
The approved on-disk file is also checked. This identifies the loaded main
executable, not its shared libraries, interpreted program, loaded scripts or
future behavior. Those limits must inform the independent service review.

Before publication it rechecks the frozen evidence, collector and boot, and
the response deadline. Raw snapshots, image results, pressure summaries and
errors are retained in `environment_worker_evidence.json`. The existing atomic,
no-overwrite writer publishes a response only after its evidence file exists.
The deadline is rechecked after serialization/hashing; a late response cannot
remain a pass. Catchable interruption leaves failed evidence, then propagates.
Failures before a valid request or an uncatchable termination may leave no
response, in which case the existing guard times out without native release.

## Validation

161 focused worker, identity, policy, observer, executor and PSI tests pass.
Tests cover file-based composition with the real `EnvironmentalReleaseGuard`:
the worker's synthetic response is read and hashed before the budget callback,
and neither component writes `go.json`. No returned dictionaries are substituted
for that file handoff. Host/process snapshots, reviewed policy and allocation
are explicitly synthetic test fixtures; this is not a native Slurm experiment.

Failure cases include unknown idle processes, missing type, observation errors,
CPU/pressure evidence failures, executable/collector errors, changed files,
deadline expiry, interrupts, incorrect request identities and overwrite refusal.
Collector controls cover missing native rows, sample errors, wrong sample index,
an advanced stream and an already released worker. A read-only local test hashes
the test process's actual loaded interpreter; a different approved file is
rejected. No ordinary host service is approved by that test.

## Remaining Production Work

The worker is not automatically started by the executor or submission script.
The final orchestration must launch it in the same allocation, retain its log
and terminal status, and clean up only its own waiting worker if preparation
fails. Validate the complete native handoff and its latency under the frozen
recipe before submitting any production identity. Existing bounded observer
captures alone do not prove it can always finish inside 20 seconds.

A real reviewed service/executable policy, prospective numerical bounds,
whole-run policy application, full-scale observer/resource validation and a
verified quiet window remain required. The latest retained observation still
shows about 93 competing cores; it is not a fresh observation from this turn.
The full 27-run order and scientific settings remain unchanged. Neither this
worker nor passing component tests establishes controlled timing or publication
readiness.

Invocation after those prerequisites, inside the intended allocation:

```sh
python -B -m benchmark_tools.threadripper_environment_worker \
  --request /absolute/request.json --request-sha256 REQUEST_SHA256 \
  --policy /absolute/policy.json --policy-sha256 POLICY_SHA256
```
