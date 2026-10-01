# Clean Test Dependency Collection

The latest scientific-code-branch [CI run 36940075541](https://github.com/JLSteenwyk/orthohmm/actions/runs/36940075541)
at `55e207e3` is terminal failure. Unlike the earlier DNS failure, installation
completed and tests reached collection. Docs succeeded. The inspected full
and Python 3.13 logs identify missing test libraries and import-time access
to Linux-only `os.sched_getaffinity`; the full job stopped with 420 collection
errors. The [measured readiness receipt](ci_collection_readiness_20261001.json)
pins both logs and the run/jobs API artifacts. Other job causes are not
inferred from those two logs. Three matrix siblings were cancelled.

Declared ten additional workflow/biology/plotting/PDF/packaging dependencies
in `tests/requirements.txt`, using versions observed in the tested development
environment. Application requirements and historical scientific runtime locks
were not changed. Both CI install steps now verify these imports explicitly.
Matrix fail-fast is disabled so one failure does not hide sibling outcomes.

The affinity observer now selects its default OS reader when called, not at
import. Without the API, real observation raises `NotImplementedError` before
reading scope files; there is no invented CPU set or successful compliance
result. Explicit synthetic readers still work. Native calibration reports its
existing capability skip rather than crashing while checking unsupported APIs.
Historical executors/receipts retain their actual prior source bindings;
this edit does not re-admit old observations or satisfy a timing handoff.

Unit, fast and coverage Make targets previously omitted the 232 top-level
test cases from CI. They now collect `tests` while excluding only the separate
integration directory; the fast target additionally excludes its existing
`slow` marker. No failing test family was removed from discovery.

## Executed Validation

- Original environment: 53 observer/calibration tests pass in 1.89 seconds.
- Fresh private Python 3.12 environment: declared requirements install,
  `pip check` passes, and all **13,681 cases collect without errors**.
- In that environment, **67 focused cases pass**, zero failures/errors/skips,
  including missing-API, injected-reader, dependency and Make-target checks.
- Structured workflow parsing and Make dry-run validate both import checks,
  matrix outcome retention and top-level discovery. Scoped whitespace passes.

The first focused command named a nonexistent third test module; pytest exited
4 without running tests. That operator selector error is retained in the receipt;
the corrected two-module command and subsequent broader command are separate.

The private environment uses newly resolved application/transitive versions,
including NumPy 2.4.6 and Leiden 0.12.0, not the frozen benchmark environment.
Collection imports current checkout sources and existing native libraries;
it is not a clean installed-package/native-runtime reproduction, a complete
new test execution or a benchmark. New pytest deprecation warnings remain
visible. No shared package upgrade, inference/scoring, contention poll, DGX,
unrelated job/service action or controlled timing run occurred.

Remote CI at the prepared revision still needs observation after the push;
the earlier failure is not reclassified. Remaining platform/runtime tests may
expose additional issues. Local collection and focused results do not establish
remote success, publication readiness, rights clearance or completion of the goal.
