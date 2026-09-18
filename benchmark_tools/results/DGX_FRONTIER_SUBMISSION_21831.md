# Native Frontier Engineering Submission

The submission SSH process exited successfully with array21831 before the
main-host clock observation2026-09-18T23:01:08Z. All transfers and the
33-file source comparison completed before submission. The launcher/audit
was committed and pushed asfd47088; the transferred inventory and identity
test were committed and pushed as9ea9afb.

Recipe SHA-256:
`1e7d4ae42107578ada71872464b6290ee6e3ae34b1a82a96af1b13e02103c7b4`.

Requested exclusive spark-7ff0,20CPU/96GiB, one-hour task limit, no requeue,
array0-2%1 and eligibility delay60seconds. This uses the frozen645-protein
fixture and the original three native commands. No DGX SSH/SCP calls are
permitted after submission completion until the local Slurm controller
shows all three tasks terminal. Local scheduler polling is permitted.

All failures and interval flags will be retained without selective reruns.
The frontier counters cannot correct native timings or relax acceptance
thresholds. This is an engineering check, not scientific timing admission.
