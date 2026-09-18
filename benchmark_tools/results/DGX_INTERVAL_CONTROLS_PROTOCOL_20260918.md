# DGX Interval CPU Controls

Freeze this engineering experiment before execution. It is not a timing
panel and does not admit historical measurements or authorize replacement
scientific runs. Use spark-7ff0 in an exclusive Slurm allocation with one
two-CPU batch task,2GiB memory and a five-minute limit, without requeue.

Run quiet then completed-burst controls. In each, a separate one-CPU native
step performs arithmetic until released. Observe eleven points on a
one-second schedule spanning ten seconds. At each point read host counters,
the native STEP aggregate cpu.stat, then host counters again. Check worker
membership before/after; keep the worker alive until the final observation.
The step aggregate includes descendants, not just the visible worker PID.

For each consecutive pair, outer host reads enclose native reads. Subtract
the native cumulative CPU delta from outer host busy CPU, retaining the
signed result and dividing by the native-read midpoint separation. Retain
all raw counters/read brackets and report outer read overhang. Require
stable host/topology/scope, no read errors, nondecreasing relevant counters,
nonoverlapping points, and intervals between0.5 and1.5seconds. Missing data
invalidate the experiment; never interpolate. Host guest time is not counted
twice and steal time is flagged separately. Outer host windows overlap, so
do not sum interval residuals to estimate total unassigned CPU.

Use the existing operational screens: residual above0.25average cores,
residual below-0.5CPU seconds, or any steal increment flags an interval.
Separately compute a whole-window screen from first/last points, without
summing intervals. These thresholds are engineering choices, not confidence
limits or rigorous bounds on interference.

The burst control starts a sibling in our batch scope after point2. It burns
0.75process-CPU seconds and must finish before point4 begins. Retain its
start/end timestamps,CPU time and cgroup. Expect quiet to pass all interval
and whole-window screens; expect the burst to flag at least one interval
while passing the whole-window average screen. Preserve unexpected results
without retuning thresholds or selectively discarding trials.

These are two fixed-order controls, not a calibration across workloads.
Observer/kernel CPU remains in residuals. Some worker startup/teardown is
outside the fixed observation window. Accounting delay, native-pipeline
overhead, CPU frequency/thermal state and memory/I/O interference are not
resolved here. No native runtime correction is permitted. Future scientific
execution still requires complete command boundaries and a frozen inclusion,
resource, run-order, missingness and repeat protocol.
