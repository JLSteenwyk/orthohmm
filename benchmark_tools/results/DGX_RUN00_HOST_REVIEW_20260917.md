# DGX Run 00 Host-Evidence Review

This is a diagnostic review, not scientific timing admission. No active run,
deployed recipe, monitor rule or result was changed. The first six array tasks
had completed while task 6 was running when this review began.

Run 00 has array identity 21656_0 and collector job ID 21657. Its retained host
summary reports `inconclusive`, despite a maximum observed persistent foreign
CPU rate of approximately 0.007534 cores, below the monitor's existing 0.25-core
detection threshold. Low observed persistent usage does not measure unmatched
or entirely unsampled processes.

Downloaded only run 00's completed host log (2,794,463 bytes) and measurement
JSON (5,511 bytes), not native inference outputs or runtime trees. This bounded
transfer took place while a later task was active and is not itself evidence of
an otherwise idle host. Raw evidence remains work-only; the committed report
records hashes and compact diagnostics without per-process command lines.

`review_dgx_host_intervals.py` verifies the raw host-log checksum against the
measurement, replays every interval through the unchanged observer, and
reconstructs its summary. All 20 snapshots and 19 intervals replay exactly:

- Seven intervals report no large persistent competitor observed.
- Twelve intervals are inconclusive.
- No snapshot errors or invalid CPU/time or changed-cgroup reasons were present.
- All 22 unmatched identity appearances/disappearances have `kworker/` names
  and root-cgroup membership in the retained samples.

Identity matching uses both PID and creation time, so PID reuse cannot silently
associate a disappeared process with a different newly created process. Names
suggest kernel worker activity but do not independently authenticate process
type. The recorded CPU contribution of unmatched workers cannot be reconstructed
from a missing endpoint. The original `inconclusive` classification is retained;
no quiet-host certification or zero-contention claim follows.

Reproduce on the downloaded evidence:

```bash
python benchmark_tools/review_dgx_host_intervals.py --directory benchmarks/work/dgx_host_review_run00_20260917 --output /tmp/dgx-run00-host-review.json
python -m pytest tests/unit/test_review_dgx_host_intervals.py tests/unit/test_command_host_monitor.py tests/unit/test_observe_host_competition.py tests/unit/test_audit_slurm_measurement.py -q
```

Evidence: `dgx_run00_host_review_20260917.json`. This closes the question of what
caused the retained inconclusive intervals for this run, not whether all timing
requirements are satisfied. Remaining admission checks include scheduler and
frozen-command contracts, native output validation, full resource replay,
host evidence for every run and repeated-run comparability. Do not restart or
discard expensive runs solely because observation was inconclusive. Retain
these limitations in the final resource analysis and explicitly distinguish
dedicated matched allocation from perfectly contention-free execution.
