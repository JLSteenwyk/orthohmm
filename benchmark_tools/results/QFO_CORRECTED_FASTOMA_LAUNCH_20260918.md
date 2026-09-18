# Corrected FastOMA Launch

Full inference job **21740** is queued with `afterok:21738` (corrected input
staging), which itself waits for successful corrected OrthoFinder admission
21731. Allocation: 180 CPUs, 720 GiB, seven days, bizon, no automatic requeue.
No corrected FastOMA inference result exists yet.

The immutable executor is
`benchmarks/work/publication_qfo_corrected_fastoma_v1` at
`7aa174c3f6bf244144c1940bfa32db2e3ca80d7b`. The retained batch script is
`qfo_corrected_fastoma_run_batch_20260918.sh`.

## Execution Contract

- Require completed 2-CPU input-staging accounting, exact staging executor,
  and the 78-proteome/984,137-protein staged manifest.
- Verify all 79 source/copy bindings, exact staged file membership and absence
  of symlinks, including the admitted corrected OrthoFinder species tree.
- Recheck pinned asset/configuration files, database, installed Java/Nextflow
  runtime trees, container-runtime binaries, Docker metadata and immutable image.
- Use `-C fastoma_corrected_execution.config`, no mutable Docker profile and no
  `-resume`. Preserve the existing scientific settings, UniProt ID conversion
  and native pair generation. Keep the prespecified collection recovery resources.
- Use a restricted environment, fixed locale, offline Nextflow 22.10.8 and
  the installed Java build. Native input, output, task work and controller
  working directories are distinct. The controller starts in a fresh directory,
  so its logs and session cache are not inherited from historical runs.
- Save the full command/environment, source identity and job metadata before
  inference. Repeat provenance/runtime checks afterward and inventory published
  native outputs, trace, logs and controller timing. Preserve failure records.
- Successful process exit is only `process_succeeded_pending_native_admission`.
  It does not admit native outputs or accuracy scores.

Output root: `benchmarks/results/qfo_corrected_fastoma_v1/`.
Staged input: `benchmarks/work/qfo_corrected_fastoma_inputs_20260918/`.
The 700-GiB task pool and 20-GiB controller allowance follow the earlier
resource configuration. Per-container resource enforcement was probed; aggregate
memory and matched timing are not established by this launch.

## Validation

All 83 focused launcher, staging, runtime and resource-probe tests passed.
They include success/failure recording, drift rejection, no implicit resume,
allocation enforcement, staged inventory checks and clean environment policy.
Bash syntax validation passed.

A real tiny Nextflow probe ran successfully with the launcher's restricted
environment. One uncached task had CPU quota/period `100000 100000` and
memory limit 268,435,456 bytes. The existing independent resource auditor passed,
and runtime identity still matched after this probe. Unmatched selectors for
biological processes are expected because the tiny workflow has only one probe
process. This did not run biological inference.

Retained evidence:

- `fastoma_clean_launch_execution_20260918.json`, SHA-256
  `7985324244236f6540fa1ac418fd0c3051cda8619c49b6f8a554041f650d6dc6`.
- `fastoma_clean_launch_probe_20260918.json`, SHA-256
  `16498766408e5747a73040be1432ea2a610e737e74fa1a830c811c9ed44a777f`.

Next: independently audit native workflow completion, expected output/ID scope,
task failures and pair semantics, then run strict pair conversion and all six
corrected QfO endpoints. This comparator uses a supplied OrthoFinder tree.
Shared-host inference is not part of the dedicated DGX timing series; installed
runtime fingerprints are not a hermetic host snapshot.
