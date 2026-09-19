# Dual Native Measurement Replay

Added `benchmark_tools/replay_dual_native_measurement.py` for the prospective
three-command diagnostic in `DUAL_BRACKET_NATIVE_PROTOCOL_20260919.md`.
No frozen collector, threshold, command, runtime or running recorder changed.

The checker requires contiguous raw point filenames beginning at zero,
recomputes both screens from raw records using the collector evaluators, and
compares them exactly with the retained report. It verifies command identity,
20 CPUs, 900-second timeout, one-second cadence, job identity, native success,
command duration, raw/embedded agreement and final native-cgroup memory.
Files are hashed before replay and rechecked afterward; changes in the raw
point inventory also fail. This is independent archive consistency replay,
not an independent implementation of the CPU arithmetic.

Focused tests include altered command/cadence/job, failure/timeout, duration,
memory scope/time/errors/counters, fabricated admission, either CPU screen,
pressure/frontier faults, missing/extra points and mid-replay evidence changes.

Run:

```sh
pytest -q tests/unit/test_replay_dual_native_measurement.py tests/unit/test_measure_native_dual_bracket_step.py tests/unit/test_probe_dual_cpu_brackets.py tests/unit/test_replay_pressure_overhead.py
```

All 63 focused tests passed in 1.09 seconds. No live diagnostic output has yet
been replayed. Scheduler success alone is
not acceptance. Separate checks remain necessary for scheduler allocation,
recipe/runtime/input hashes, authorized task selection, complete native outputs
and canonical output equivalence with the retained earlier panel. Neither
CPU screen proves foreign-workload absence. Cgroup peak memory includes
wrappers/cache and is not maximum process RSS. Observation windows are not
exact command boundaries. Comparative timing admission remains false.
