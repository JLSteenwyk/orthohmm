# Prospective DGX Counter Availability

## Evidence

Executed a read-only three-second probe on spark-7ff0 over the authorized
Ethernet SSH connection. Kernel6.11.0-1016-nvidia exposes the aggregate CPU
row, online CPU list, boot identity, CPU/memory/I/O pressure, and observer
cgroup cpu.stat, memory.current, memory.peak and memory.events. Both snapshots
read all requested fields without errors. Raw text and monotonic read brackets
are preserved in dgx_host_counter_probe_20260918.json. Local replay exactly
reproduces its derived summary. The source SHA-256 is
3d3c193305085e17aa10e05be842d814f66d020d7dde285c8d9288ead2ce3f3c.

The observed midpoint interval was3.000967137 seconds, with0.01 accounted host
busy CPU seconds. This is a short availability diagnostic, not an exclusivity
certificate or a prediction of load during future runs. It includes SSH and
the probe itself. The scope is the SSH session's cgroup, not a Slurm job or
native inference subtree. No foreign-CPU subtraction or timing correction is
made. No unrelated process was stopped and no scientific run was repeated.

## Accounting Semantics

The [kernel proc documentation](https://docs.kernel.org/filesystems/proc.html)
defines aggregate CPU fields in USER_HZ units and warns that iowait can
decrease. The probe obtains SC_CLK_TCK rather than assuming a tick rate,
retains iowait decreases explicitly, and rejects decreases in other fields.
Busy CPU sums user,nice,system,irq,softirq; guest fields are not added again,
and steal remains separately reported. Non-atomic midpoint rates are not
rigorous bounds on workload interference.

[Cgroup v2 documentation](https://docs.kernel.org/admin-guide/cgroup-v2.html)
describes cumulative CPU accounting and cgroup memory controls.
[PSI documentation](https://docs.kernel.org/accounting/psi.html) describes
stall pressure. These counters supplement process snapshots; they do not
identify arbitrary foreign workloads or prove absence of interference.

## Required Validation Before Timing

1. Capture these fields prospectively inside a scheduled job, verifying the
   intended native/observer cgroup boundaries and node identity. An SSH-session
   probe is not that validation.
2. Use a separately frozen engineering test with known short-lived external
   CPU work that exits between process scans, plus a quiet control and sustained
   load. Verify what host and cgroup counters actually retain. Retain failures.
3. Measure observer overhead and accounting discrepancies under compute-heavy
   and process-heavy workloads. Do not treat host minus cgroup CPU as a proven
   foreign-load upper bound: read intervals and accounting scopes differ.
4. Define new scientific inclusion rules, cadence, memory scope and repeat
   policy before examining any replacement benchmark times. An observer change
   cannot retroactively supply the original panel's missing observations.

The existing27 runs remain descriptive-only. A full automatic rerun is not
authorized by this probe. Eighteen focused tests cover counter arithmetic,
guest/steal handling, iowait decreases, changed identity/topology, malformed
records, counter rollback and a local read-only smoke. These tests and the
DGX probe establish availability and arithmetic only, not complete observer
validation or publication-ready controlled timing.
