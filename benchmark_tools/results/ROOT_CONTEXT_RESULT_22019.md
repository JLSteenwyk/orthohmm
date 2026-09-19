# Root Context Panel 22019: Incomplete

## Outcome

Slurm job 22019 failed with exit `1:0` after 1 minute 45 seconds. It was
submitted once after a recorded empty DGX queue, with exclusive 20 CPUs,
96 GiB and a 15-minute limit. No SSH inspection or transfer occurred between
submission and the terminal scheduler observation. Local controller capture
completed with seven polls and zero observation errors.

The whole-panel audit reproduces source/runtime bindings, scheduler limits,
all 12 outcome entries and cumulative checkpoints. It independently validates
the first three native workloads and their raw observations. It does not
validate the failed trial or count unrun conditions as successes.

| Index | Condition | Outcome | All/common intervals | Original/narrow flags |
| --- | --- | --- | --- | --- |
| 0 | Idle | Validated | 21 / 19 | 0 / 0 |
| 1 | Steady | Validated | 21 / 19 | 0 / 0 |
| 2 | Churn | Validated | 21 / 19 | 0 / 0 |
| 3 | User-contended | Execution failed | Not admitted | Not admitted |
| 4-11 | Remaining fixed conditions | Unrun after failure | Missing | Missing |

Aggregate user-slice CPU over each valid enclosing observation span was
0.033673, 0 and 0 seconds respectively. These are not isolated task CPU
measurements. No user-contended positive control completed, so response
sensitivity remains untested. All signed residuals and descriptive within-block
differences are retained, with missing comparisons left null. Dependent
intervals are not treated as independent replicates.

## Failure Evidence

Trial 3's service log reports `Failed to connect to bus: No such file or
directory`. The coordinator did not receive `competitor_ready.json`, did not
release the shared workload start, and recorded a timeout. Native workers
failed while waiting; the panel stopped all remaining launches.

Post-termination journal retrieval shows the user manager starting at
22:42:38 UTC and shutting down at 22:42:48 UTC, ten seconds after submission.
A subsequent `loginctl show-user 1000` inspection reported `Linger=no`.
This explains why the earlier preflight could observe the manager but the
later owned-service command could not connect. This is a session-lifetime
failure, not evidence that the intended CPU load was missed by the probe.
No unrelated services or persistent login configuration were changed.

Before a new experiment, document and validate a bounded way to maintain
the required user-manager session for the entire panel. Preserve this failed
panel; any subsequent full run must be identified separately, not spliced
into its missing repetitions. Do not silently enable persistent user
services or reinterpret failed controls as low responses.

## Retained Evidence

- `root_context_archive_22019.tar.gz`: all 805 retrieved files, including
  raw failed evidence, source recipe tree and job log; SHA-256
  `65c9f2a65c1dadf5f4738c4873b1e7a0c5f62b9768b027c8bfcd74ea76996147`.
- `root_context_scheduler_22019.tar.gz`: complete local scheduler capture;
  SHA-256 `0c52e49474d73a8e45c3524a0f54348910a520a7b7850d7ab6e4482c0094a130`.
- `root_context_audit_22019_20260919.json.gz`: full audit; SHA-256
  `06bdeb7334454c4e616fbf783334d33bd5c72b3b3e73aafd7e452f87d24f639f`.
- `root_context_description_22019_20260919.json`: all/common interval
  distributions, flags, root membership changes and fixed-block differences;
  SHA-256 `6badc7db1cf47120bda31258bc8a4cdc82a101ded0553e7eb0c2dbc88dfa7466`.
- `root_context_user_manager_journal_22019.json`: bounded historical journal
  query and raw response; SHA-256
  `a82fb267fd6e80880478f0856a7347da5b71c67705fb9edfd3f234cec1d367ff`.

The archive can be extracted into a fresh directory and replayed using
`benchmark_tools.audit_root_context_controls` with job 22019, the retained
recipe, its pinned SHA from `ROOT_CONTEXT_DEPLOYMENT_20260919.md`, and the
terminal scheduler record. Description uses `describe_root_context_controls`
with the audit SHA above. No native overhead, attribution specificity,
non-CPU isolation or scientific timing eligibility is established.
