# Native QfO Measurement Failure: Original 22437

Original index7/P0C0R1 completes its native Slurm step with exit0:0 in
12:55:21. The enclosing allocation fails1:0 after13:04:30, and original
reviewer22440 fails1:0 after34seconds. These are distinct retained outcomes,
not a clean-success timing. No inference or review is restarted.

## Cause And Scope

Original `verification.json` and `measurement/report_finalization.json`
report `Missing or irregular observation interval`. Finalization fails
before creating `lineage_report.json` and `root_context_report.json`.
The original reviewer completes its runtime review, then fails when the
unchanged replay requires the absent lineage report. Original receipts,
raw observations and scheduler outcomes remain untouched.

The [new read-only cadence audit](native_factorial_cadence_failure_22437_20261006.json)
uses [this source](../audit_native_factorial_cadence_failure.py) and the
unchanged interval API. It streams all46,511 numbered points once,
9,991,918,809 raw bytes, and retains every interval rejection:

| Left Point | Right Point | Observed Period (s) | Unchanged Criterion |
| --- | --- | ---: | --- |
| 23210 | 23211 | 1.6594024915 | Above1.5s maximum |
| 23211 | 23212 | 0.3398767975 | Below0.5s minimum |

No other interval-API rejection occurs among46,510 intervals. This is
not a complete lineage, affinity or resource replay. Individual hierarchy/
interval point validation and unchanged pair rejection checks are not a
substitute for all scientific measurement requirements. The ordered
path/bytes/SHA256 inventory digest is
`0c81b6de2bd8f6a75b317215f4f16124db0a0edaade929ccfe36ba1fad7ab0ab`.
Only boundary/failing-point references are embedded; the full raw data remain
in their original persistent location, not committed as a large dataset.

The retained pressure review separately fails evidence cadence on one
1.659820477s period. Its diagnostic pressure exceedances are not converted
into exclusion grounds for ordinary contention. The pressure and primary
intervals have different read brackets; their periods are not interchangeable.
The cause of the delayed observation is not established. No interpolation,
threshold relaxation, time correction or contention-only retry is performed.

All920 prospectively frozen helper records still match. The audit is
11,067bytes/SHA256
`78a13590e1d33e89c91a80c88ad41e561d1475f2c1ce1cbb1663df5a99949e99`.
It deliberately keeps full-resource replay, scientific timing admission,
scientific-output validation and successor authorization false; primary
resources are null, not zero or invented.

## Scientific Completion And Next Action

Original native receipt reports completed-pending-output-review, native
metrics report complete, and `done.json` records exit0/no timeout. Metrics
report397,041 final groups and5,113,820 resolved ortholog pairs. Those are
completion records, not independently admitted predictions or accuracy.
No new QfO endpoint, bootstrap interval or historical output equivalence is
asserted. Existing admitted first-cell scores remain unchanged.

Next perform an explicit, separately recorded failure disposition and
scientific semantic-output review. Reuse the original successful runtime
review only after checking its bindings. Keep invalid-cadence resource
evidence excluded from eligible timing comparisons. A successor requires
an explicit reviewed disposition plus its unchanged fresh safe-capacity,
input/runtime and sequential-order gates; this diagnosis grants none.
Retain the failed original reviewer instead of silently repeating it against
fabricated reports. Scoring must separately bind recovered scientific
outputs and disclose the failed measurement route.

## Validation

[25 joined tests](native_factorial_cadence_failure_tests_20261006_v2.xml)
pass in0.63s, zero failures/errors/skips, including10 new cadence-audit cases.
They cover exact original boundaries, multiple slow/catch-up rejections,
non-cadence errors, invalid points, missing observations and non-admission.
The [initial collection error](native_factorial_cadence_failure_tests_20261006.xml)
is retained: the new test used an unqualified sibling import. Correct only
that new test import to the repository's `tests.unit` pattern; production
and frozen measurement code do not change.

Bounded independent readback checks the three failing-interval raw points,
recomputes both periods with exact rational arithmetic, verifies12 direct
evidence bindings and reads the25-test XML. It also confirms the new helper
is outside the920 frozen sources and all non-admission flags remain false.
The9.99GB census is not repeated for this readback.

An unscoped whitespace inspection also found pre-existing whitespace in
unrelated generated samples. None of those files are changed or committed.
The focused source checks pass. The audit is diagnosis, not publication
readiness, independent scientific validation or an isolated-performance claim.
