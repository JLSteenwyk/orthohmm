# Corrected Proteinortho Submitted

Job **21708** was submitted from detached executor
`benchmarks/work/publication_qfo_corrected_proteinortho_v1`, revision
`989f5d9185018f2f467a25842a665dad08b92910`. Scheduler inspection confirmed
RUNNING on bizon with 32 CPUs, 192 GiB, 72-hour limit, no requeue and zero
restarts. No original-release output directory was reused.

The execution record under
`benchmarks/results/qfo_corrected_proteinortho_v1/execution.json` confirms
78 copied input files, running status, and these exact provenance bindings:

- Command plan: `65f9baedeb6f6cbb66198e2a6a24ae8531b48eea9a04809c2c2443d570d4f93a`.
- Runtime snapshot: `e0fc669cbbd805f9ec768ccd9088977609be3b4a88d9ff95bfb74f5b73f33ae5`.

The native log independently reports Proteinortho 6.3.6, DIAMOND 2.1.12,
32 CPU threads, and successful initial input checks. Preflight runtime
checks passed both before staging and from the actual native working
directory. The initial cache-directory `find` warning also occurred in
the historical run; it is preserved, not removed from the log. Successful
startup does not establish successful completion.

Native outputs, conversion and accuracy remain unadmitted. This run is
shared-host accuracy work, not matched timing evidence. The ongoing
OrthoHMM/OrthoFinder corrected primary array, original-release factorial
jobs and dedicated DGX timing array are unchanged. Remaining corrected
comparators and factorial cells remain required.

Full local unit suite after launcher implementation: **2,612 passed in
56.07 seconds**. This verifies the software test suite, not unfinished
native benchmark outputs or publication readiness.
