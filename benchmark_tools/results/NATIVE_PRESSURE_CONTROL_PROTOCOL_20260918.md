# Prospective Native CPU Pressure Controls

Freeze before collecting control outcomes. Purpose: check whether the native
step's CPU PSI responds to known overlapping CPU demand on the same allocated
core. This is not a scientific timing experiment or calibration of arbitrary
interference. No OrthoHMM or comparator inference runs in this experiment.

Use a two-CPU/256MiB Slurm job on idle spark-7ff0, no requeue, five-minute
limit. The batch observer and native step remain distinct. The worker pins
itself to one CPU in its scheduler-provided affinity; the observer pins to a
different permitted CPU. A batch-step competitor, when used, explicitly pins
to the worker's CPU within the same job allocation. No host services or
unrelated jobs are stopped or moved.

Nine sequential trials, three of each mode. Fixed block orders:

1. quiet, native-only, contended
2. native-only, contended, quiet
3. contended, quiet, native-only

Quiet: native worker sleeps 1.5 seconds. Native-only: worker executes the
existing `probe_dgx_step_separation.burn` routine targeting 0.75 process-CPU
seconds. Contended: the same native work overlaps a batch-step subprocess
executing the same burn on the same CPU. Retain affinity, memberships,
process-CPU durations and monotonic work boundaries for both participants.
After native completion, wait 0.2 seconds before the final pressure read in
every mode; retain the worker until the after-observation is complete.

Save before/after native PSI with host brackets for all trials. Use the
already tested scope/identity and monotonic-counter validation without
changing thresholds in the older timing panels. Independently replay the
new report and source identities after collection.

## Fixed Engineering Checks

- Native and competitor work must each report 0.75-0.90 process-CPU seconds
  when present. Require successful worker exit and unchanged intended
  affinities and step identities.
- In each contended trial, recorded native/competitor work windows must
  overlap for at least 0.25 wall seconds. Invalid injection is not a negative
  sensitivity result and must remain explicitly failed.
- Within each block, contended native CPU `some` stall total minus
  native-only total must be at least 100,000 microseconds. Require this
  response in all three blocks for the predefined response check to pass.
- Report quiet totals and all native CPU `full`, memory and I/O totals as
  diagnostics. Do not demand zero baseline pressure or infer memory/I/O
  detection sensitivity from a CPU-only workload.

Keep every trial and failure. Do not extend individual burns, change core
selection, choose favorable blocks, rerun failures selectively, or derive an
inference exclusion threshold from these outcomes. A failed expectation
requires diagnosis and a separately documented future experiment.

Even a passing result validates only this injected contention response under
this configuration. It does not measure observer overhead, explain residuals
in previous panels, establish host isolation, or authorize any of the 27
scientific scaling runs. Overhead and scientific inclusion remain separate
prospective requirements.
