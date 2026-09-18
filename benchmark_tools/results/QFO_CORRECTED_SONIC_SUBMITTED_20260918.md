# Corrected SonicParanoid Submission

Submitted job **21710** from detached executor
`benchmarks/work/publication_qfo_corrected_sonic_v1`, revision
`d3f6dc0a7931579af779470b7f982e34a9d09b2b`.
Scheduler inspection confirmed RUNNING on bizon, 32 CPUs, 192 GiB,
72-hour limit, no requeue and zero restarts. The job initially performs
full preflight; RUNNING alone does not establish native inference startup
or successful completion. The standalone real preflight passed before
submission, as recorded in QFO_CORRECTED_SONIC_LAUNCHER_20260918.md.

Frozen plan:
`188efb1f18e8f0bd90967670d18015ad500a9c47480333b0b380dff7c10c4c59`.
Frozen resolver snapshot:
`b4ed57defaba417f6627b671df958a2d533f4438bb8ff35a0b95cb820fd61bc9`.
Output root: `benchmarks/results/qfo_corrected_sonic_v1`.
Scheduler log: `benchmarks/work/qfo_corrected_sonic_21710.log`.

Next inspect the scheduled preflight and native log, then terminal native
output completeness and species-pair provenance before any conversion or
scoring. No SonicParanoid accuracy result is yet admitted. This is
shared-host accuracy work, not dedicated timing. Existing Proteinortho,
OrthoHMM and original-release factorial jobs remain untouched.

DGX task 21656_13 independently completed with exit 0:0 in 35:13; task
21656_14 was observed running. These are scheduler states, not completed
native/resource/quiet-host admission for the timing series. Remaining
corrected comparators, all corrected factorial cells and publication
deliverables remain in scope.
