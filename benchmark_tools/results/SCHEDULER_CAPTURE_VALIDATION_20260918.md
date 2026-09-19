# Automatic Detailed Scheduler Capture

The failed evidence collection for array 21838 is not repaired retroactively.
New helper `benchmark_tools/capture_array_scheduler.py` at committed revision
`25ff71258b06e30974ae62512ac1d0066a20031c` automates future controller-side
collection before terminal records expire. It never queries the compute host
over SSH, substitutes `sacct` summaries, or admits scientific timing.

## Tests And Live Validation

Twenty collector tests and 38 existing provenance tests pass (58 total).
They cover sequential completion followed by controller expiry, missing
records, transient timeouts, malformed records, accounting-only text,
compressed pending arrays, failed tasks, immutable first terminal records,
and refusal to overwrite an existing capture directory. The expiry scenario
is simulated; the live test does not change Slurm retention configuration.

Submitted local smoke array 21863 on `bizon`, not the DGX. Both tasks only
executed `sleep 3`, each with one CPU and 64 MiB requested memory. Both
completed 0:0 in three seconds. Collector exit 0 retained both detailed
terminal records at polls 20 and 50, across 51 polls and zero observation
errors. No unrelated job or service was changed.

The [capture receipt](scheduler_capture_smoke_21863/capture.json) has SHA-256
`9d654e40cbece031ebb1cba88ca4505c65345104871b48d4a9f9e5629a83572b`.
Its adjacent `scheduler_0.txt` and `scheduler_1.txt` are unchanged copies of
the raw terminal records; receipt hashes match both. All raw poll responses
remain in `benchmarks/work/scheduler_capture_smoke_21863/`.

Commands used at repository root:

```bash
sbatch --parsable --hold --job-name=scheduler_capture_smoke \
  --partition=gpu --nodelist=bizon --array=0-1%1 --cpus-per-task=1 \
  --mem=64M --time=00:02:00 --no-requeue \
  --output=benchmarks/work/scheduler_capture_smoke_%A_%a.log --wrap='sleep 3'
/home/bizon/anaconda3/bin/python benchmark_tools/capture_array_scheduler.py \
  --array-id 21863 --tasks 2 --output benchmarks/work/scheduler_capture_smoke_21863 \
  --interval 1 --max-seconds 180
```

While the collector ran in its own terminal session, `scontrol release 21863`
released the held array. A new smoke test must use its newly assigned ID and
a fresh output directory, not reuse 21863.

## Future Timing Workflow

Submit a future frozen array held, start this collector on the controller
host, verify it has written a successful initial observation, then release
the array. Use a collection deadline covering the entire serial array, not
only one task. Default polling is five seconds; choose an interval well
below the controller retention period. Retain the collector process handle
and check for exit or errors while jobs run. Collection can still miss
records during extended outages; it reports an incomplete result rather
than reconstructing missing fields.

After the array is terminal, preserve `capture.json`, all raw polls and
`scheduler_N.txt` files with the measurement archive. The independent audit
must still validate exact allocations, commands, restarts, source identity,
outputs and resource evidence. A complete capture can include failed jobs;
it is not a successful benchmark panel.

Transient cgroup handling, non-CPU isolation, and scientific timing inclusion
remain unresolved. This collector alone does not justify another scaling
panel, validate monitoring overhead, or admit any of the original 27 runs.
