# Prospective Dual-Bracket CPU Controls

## Scope

Freeze before collecting new outcomes. This experiment tests the earlier
host bracket as a CPU-accounting diagnostic alongside the unchanged wider
bracket and native PSI. It does not replace or repair panel 21889, admit
scientific timing, or calibrate arbitrary interference. Historical flags
and failed tasks remain unchanged.

Use `probe_dual_cpu_brackets.read_point` and `compare`, which delegate
collection and threshold arithmetic to the existing frozen helpers. Retain
full raw points, including both host endpoints, frontier counters and native
pressure. No hardware-counter interpolation, residual clipping or wall-time
correction is allowed. Preserve signed discrepancies.

## Fixed Design

Use a dedicated DGX Slurm allocation: two CPUs, 256 MiB, exclusive node,
five-minute limit, no requeue. Verify no other scheduled workloads before
submission. Do not stop, move or modify unrelated services or jobs.

Run nine sequential trials with the same workload definitions and balanced
orders used in the prior native-pressure controls:

1. quiet, native-only, contended
2. native-only, contended, quiet
3. contended, quiet, native-only

The sleeping/CPU-burning native worker occupies a distinct one-CPU Slurm
step and remains alive through final measurement. Pin the observer to a
different permitted CPU. Quiet sleeps 1.5 seconds. Native-only and contended
use the existing 0.75-process-CPU-second burn; accept doses only within
0.75-0.90 seconds. The contended mode adds the same burn from a batch-step
subprocess pinned to the native worker's CPU. Require recorded overlap of
at least 0.25 seconds, correct memberships and affinities, and zero exits.
Wait 0.2 seconds after both workloads finish before the final read.

Take one full dual-bracket point before release and one after completion.
These are whole-control windows, not periodic samples: explicitly disable
the periodic 0.5-to-1.5-second gap restriction, retaining strict chronological
ordering and all other validations. Retain native and competing work times
and confirm both are enclosed by the measurement window.

## Prespecified Checks

- Validate injection identity, dose, overlap and observation enclosure before
  interpreting a response. Invalid injection is failed evidence, not a
  negative sensitivity result.
- Report all wider and narrower residuals and flags using unchanged limits:
  positive residual above 0.25 average cores, negative discrepancy below
  -0.5 CPU-seconds, or nonzero host steal.
- Require narrower screens to pass in all three quiet and all three
  native-only trials, and to flag positive excess CPU in all three valid
  contended trials. Retain wider-screen classifications regardless.
- Require native CPU PSI `some` in each contended trial to exceed its
  block's native-only total by at least 100,000 microseconds, as in the
  previous pressure protocol. Report full, memory and I/O totals without
  inferring their sensitivity from this CPU-only experiment.
- Keep all nine trials and any failure. No selective repeats, dose changes,
  alternative CPU selection after outcomes, or threshold relaxation.

## Evidence and Limits

Before dispatch, pin the protocol, runner, reader and transitive helper
source hashes in a fresh recipe; retain exact Python/runtime and scheduler
identity. Independently replay raw points, injection validation and summary
after the entire job terminates. A missing or failed trial prevents a
complete passing control result.

Even success establishes only single-core workload response under this
configuration. It does not resolve full-node/native-tool bursts, cgroup
identity churn, memory or I/O interference, observer overhead, or matched
scaling eligibility. Those require separate prospective evaluation, with
the complete overhead panel repeated rather than only selected pairs.
