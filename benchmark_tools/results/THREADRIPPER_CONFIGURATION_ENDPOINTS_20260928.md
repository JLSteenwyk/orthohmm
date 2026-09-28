# Explicit Configuration Checks At Post-Run Review

The environmental worker already checks the declared service/configuration
files and includes them in its preflight evidence. Previously, post-run review
rechecked these files only when that evidence list carried the references. It
did not itself require the policy's configuration inventory or directly bind
its entries. This is a contract-hardening change, not evidence that a real
production timing was previously admitted with changed configuration.

The bound stream audit now requires a nonempty `configuration_files` list and
adds every declared reference to the evidence checked before and after review.
The result explicitly records `configuration_endpoint_hashes_verified` only
after those checks pass. Missing declarations, empty/non-list inventories and
changed bytes fail rather than producing an accepted review artifact.

Six selected tests failed before the change, including the new positive-result
field and five missing/changed-evidence rejection cases. After the change, all
168 focused stream, pressure, process-policy, environmental-worker and executor
tests pass. File-backed mutation tests cover changes before review and during
process-stream evaluation without relying on incidental preflight references.

## Boundary

These checks compare endpoint bytes to reviewed hashes. They do not observe a
transient configuration change restored before review, attest runtime settings
loaded from other sources, prove inventory completeness or check loaded shared
libraries. Whole-run executable/configuration stability remains unresolved.
The updated helper must be included in the next execution-source freeze and
native/full-scale validation; test success does not authorize timing submission.

## Fresh Host Observation

Slurm was empty, but a fresh three-second typed process capture measured
92.4384 competing CPU-core equivalents, with six sampling errors and seven
unmatched processes retained. BALiPhy, matched-interval and Neocallimastix guide
services each accounted for approximately 16 cores. Several other scientific
services each used approximately eight cores. The
[retained summary](threadripper_timing_recheck_after_wgd_20260928.json) binds
the raw capture and lists the eight largest observed cgroups.

This is a dated persistent-process observation, not a complete measurement of
host demand or a verified quiet window. No production benchmark was launched,
unrelated job/service modified or DGX accessed. Existing quiet-window
coordination remains necessary; no timing result was admitted.
