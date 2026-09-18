# Corrected FastOMA Native Admission

Independent validation job **21741** is queued with `afterany:21740`, the
corrected FastOMA inference job. It reserves 2 CPUs and 64 GiB for four hours
on bizon, without automatic requeue. It cannot admit failed inference: the
validator first requires terminal `COMPLETED`, exit `0:0` accounting for the
180-CPU/720-GiB inference allocation.

The immutable executor is
`benchmarks/work/publication_qfo_corrected_fastoma_admission_v1` at
`857bb5e0d6adad9960ea58d98a3037097ee40d5f`. The retained batch script is
`qfo_corrected_fastoma_admit_batch_20260918.sh`.

## Checks

- Exact frozen inference source and freshly repeated staging/runtime preflight;
  command, environment, working directory, job identity and finite timestamps.
- Complete published output inventory matching the successful execution record.
- Strict native task coverage and container-wrapper checks. Failed/retried or
  cached task chains require separate review rather than silent acceptance.
- Published XML, root/marker tables, checked tree, report and native pair file
  matching their producing task artifacts. Pair extraction must consume the
  same XML that was published.
- Each OMAmer query matching its staged corrected proteome; the resulting map
  matching both the published map and the root-inference input map.
- Exact 78-species/984,137-protein input universe. The checked species tree
  must preserve supplied rooted clades; internal-label renaming and child order
  are allowed, but altered topology, species scope and nonfinite branches fail.
- Native XML identifier/taxonomy consistency, exact root-table membership and
  cross-species pair membership within root HOGs. Unassigned protein counts
  remain explicit rather than being treated as complete coverage.
- Rehashed input/source/output/task records and unchanged published membership
  before a success report is written.

The prospective report is
`benchmarks/work/qfo_corrected_fastoma_admission_20260918.json` with status
`corrected_fastoma_native_evidence_admitted`, `accuracy_evaluated: false` and
`publication_ready: false`. Existing reports are never overwritten.

## Validation And Limits

All 113 focused native-admission, task, XML and launcher tests pass. Bash syntax
validation passed. Tests cover scheduler/source/environment drift, invalid
timestamps, published-file substitution, query/map substitution, tree scope,
task coverage and XML/pair semantics. Components have also been exercised on
the real historical XML/pairs and tiny fresh container probe in the preceding
audits. **Actual corrected native admission has not run yet.**

This is an integrity gate, not proof of biological accuracy. It does not claim
that root-HOG co-membership proves orthology, that all scientific arguments were
independently reimplemented, or that shared-host resources are matched timing.
Successful admission must be followed by strict distinct-pair conversion,
reference mapping, all six QfO endpoints and independent score validation.
The supplied tree remains explicitly attributed to corrected OrthoFinder.
