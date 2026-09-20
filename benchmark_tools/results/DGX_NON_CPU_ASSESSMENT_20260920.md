# DGX Non-CPU Environmental Assessment

## Retained Observations

The pinned, validated 22022 overhead audit already retains native-step CPU,
memory and I/O pressure over each measurement observation window, plus final
memory counters. `describe_root_context_non_cpu.py` extracts all 18 outcomes
without repeating inference, changing flags or defining new inclusion gates.
The machine-readable result is `root_context_non_cpu_22022_20260920.json`.
Its SHA-256 is
`fc996509adda1fdab9e339970dfce5cc04bb397658ee8657f8d5811232d2e2a9`.

All 18 runs observed positive I/O pressure totals. Native-step `some` I/O
stall totals range from 19,507 to 1,143,527 microseconds; `full` totals range
from 14,776 to 1,137,778 microseconds. Twelve runs observed positive memory
pressure: the six satellite_v2 and six full OrthoFinder runs. Their `some`
totals are 1-392 microseconds and `full` totals are 1-122 microseconds. All
six high-sensitivity runs observed zero memory pressure in this window.
All recorded memory-event counters, including limit and OOM events, are zero.
Native-step peak charged memory ranges from 2,220,634,112 to 7,195,627,520 bytes.
All per-run values and all 344 original/60 narrow flags remain in the result.

These are observation-window counters, not exact native-command-only stall
times. CPU, memory and I/O pressure are not added together. `some` and `full`
are overlapping categories, not additive quantities. The kernel defines
pressure as task stalls, not a measure of the identity of competing workloads:
[kernel PSI documentation](https://docs.kernel.org/accounting/psi.html).
Charged cgroup memory includes file cache and kernel allocations and is not
maximum-process RSS:
[cgroup-v2 documentation](https://docs.kernel.org/admin-guide/cgroup-v2.html).

## Current Environment Check

A read-only September 20 DGX service query confirms
`samwise-daemon-samwise.service` is enabled, with `Restart=always`, a five-second
restart delay and an auto-restart state following an exit-code failure. Its
configured working directory is `/home/jlsteenwyk/Desktop/Samwise`; the retained
22022 journal records missing-working-directory failures. The current query
is retained locally as `benchmarks/work/dgx_samwise_status_20260920.txt`.
Its SHA-256 is
`1f4ebbbdf833b709fb9ec0e38d2b04a986e6f5729ce407c88fbf08b640de4a3b`.
New login sessions can create a new user manager, so a current restart counter
of zero does not contradict the 2,297 restarts in the archived job window.

A spot `nvidia-smi` query reported NVIDIA GB10 utilization 0%, temperature 41 C
and power 3.53 W; memory-used/total fields were unavailable (`[N/A]`). This is
a transient observation, not evidence about GPU activity during 22022 or a
future timing run. No unavailable field is interpreted as zero memory use.

Approval was requested to temporarily stop and runtime-mask only the failing
user service during timing, restoring its prior configuration afterward. No
approval has been received and no service has been changed as of this entry.

## Prospective Timing Decision

The completed overhead experiment supports its bounded incremental collector
comparison. It does not establish device-level I/O attribution, GPU/memory-
bandwidth isolation, thermal behavior or a quiet future execution window.
Historical resource measurements remain descriptive and unchanged.

Do not make zero native PSI an inclusion requirement: native algorithms can
cause their own memory/I/O stalls, and this would select outcomes by method.
Do not subtract native from host PSI or identify outside activity from a
signed CPU residual. Neither operation supplies defensible causal attribution.

Before launching a new matched scaling panel, establish a recorded service
and workload policy for the DGX, including the response to an unauthorized
background job or service change. Retain pre/post hardware and service state
and whole-run workload observations, report measurement limits explicitly,
and freeze handling of failures/contamination before inspecting method-speed
comparisons. The current data do not authorize retroactive admission or
alter the frozen scientific endpoints. No new scaling inference was launched.

## Reproduction

```bash
python -B -m benchmark_tools.describe_root_context_non_cpu \
  --audit benchmark_tools/results/root_context_overhead_audit_22022_20260920.json.gz \
  --output benchmarks/work/root_context_non_cpu_reproduction.json
```

The reader pins the full compressed audit SHA-256, checks source bytes before
and after reporting, rejects incomplete outcomes and invalid counters, and
retains positive stalls/events without certifying isolation. Regression tests
reproduce the description from a relocated copy of the actual audit. This is
post-outcome environmental description, not a calibrated interference detector.
