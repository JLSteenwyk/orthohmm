# Bracketed CPU Control Results

Job21805 completed0:0 in8seconds with zero restarts, exclusive spark-7ff0
allocation, two CPUs/task and2GiB requested. Exclusivity allocated20CPUs;
each native control requested one CPU. Protocol and runner were committed
and pushed as2dac972 before submission. The launcher checked all four Python
source hashes and the protocol hash before execution.

| Control | Signed unassigned CPU seconds | Unassigned average cores | Operational result |
| --- | ---: | ---: | --- |
| Quiet | 0.019355 | 0.009674 | Passed |
| Completed one-CPU-second sibling burst | 1.049399 | 0.524520 | Excess unassigned CPU |
| Sustained three-CPU-second sibling load | 3.059178 | 1.529158 | Excess unassigned CPU |

All three frozen expectations held. The outer host snapshots fully enclosed
the native snapshots; the latter enclosed the two-process-CPU-second native
work. Both native and sibling durations met their frozen ranges. Native and
observer steps differed, while each sibling remained in the batch scope.
All requested counter reads succeeded. Both positive loads had exited before
the host post-snapshot. No threshold was changed and no trial was retried.

The sustained control intentionally includes CPU executed after native work
if the sibling outlasts it. Its residual is therefore not a claim about the
amount of contention during native inference. The residual also includes
observer, kernel and launch/exit activity; no foreign process is identified
by the arithmetic. These results demonstrate coarse detection under this
fixture, not a calibrated interference bound or scientific timing admission.

## Evidence

`dgx_bracketed_controls_21805.json` retains all raw snapshots, commands,
durations, source hashes, operational results and control expectations.
`dgx_bracketed_scheduler_21805.txt` retains terminal controller text verbatim,
including trailing whitespace. Complete handshakes and logs are retained in
`benchmarks/work/dgx_bracket_probe_21805/`; remote recipe:
`/tmp/orthohmm-bracket-probe.bnlIY5VJ/`.

Local replay of the three operational results exactly matches retained
values. Regression tests also check committed source identities, actual
scheduler completion, rejected duration/scope/boundary mutations, and the
unchanged false admission flags. The previous native smokes and historical
27-run panel are not retroactively upgraded.

## Remaining Work

The corrected handshake is now exercised on known CPU controls. Interval-level
counter observation is still needed to avoid concealing concentrated bursts
in long-run averages. General observer overhead, non-CPU contention, memory
accounting and prospective scaling inclusion/repeat rules also remain to be
addressed. This experiment submits no scientific scaling runs and makes no
OrthoHMM/OrthoFinder speed or memory comparison.
