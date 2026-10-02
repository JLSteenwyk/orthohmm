# Test Security And Observed Numerical Replay Difference

Previous goal turn made progress at pushed `1b15b57d`. Reread the objective
and ledger; observe its actual [CI run 36949483608](https://github.com/JLSteenwyk/orthohmm/actions/runs/36949483608)
without resubmitting it. Linux CPU-wheel/docs succeeded; the inspected Python
3.12 job is terminal failure while siblings were still live at that snapshot.
Neither installation success nor one inspected log clears the full matrix.
At the later 2 October 01:23:23 UTC checkpoint all four fast jobs have failed;
full-suite job 110658992661 remains live. Only the Python 3.12 log is inspected
here. Preserve that full-suite handle and do not infer sibling failure causes.

## Actual Replay Difference

The [retained diagnostic receipt](ci_exact_replay_difference_20261001.json)
pins the downloaded Python 3.12 log: **13,395 passed, 189 failed, 38 errors,
96 skipped in 421.20 seconds**. All ten actual frozen-source archive tests
now pass, confirming history availability for their Git operations. All seven
count-level matched-graph reproduction tests pass under their existing 1e-12
tolerance. These counts are not like-for-like with an earlier Python 3.13 job.

The renderer reports exactly one unequal leaf in its checked contrasts and
bootstrap sections: `uneven_taxa.recall.marginal_95_percent_ci[0]` is retained
as `-0.7992007992008033` and recomputed as `-0.7992007992008032`. The absolute
difference `1.1102230246251565e-16` is exactly one ULP. No checked F1 or bootstrap
metadata differs. This is an exploratory recall interval, not the plotted F1
endpoint. No CPU-kernel or fused-operation cause is isolated. Preserve the
strict guard and historical result; a follow-up must explicitly reconcile
portable numerical admission with the count-level reproduction contract,
not silently rewrite scores or weaken checks merely to obtain green CI.

## Test-Pin Remediation

A fresh read-only API snapshot has **55 open alerts** on five manifests.
Fifteen (ten high, five medium) concern current test requirements; forty concern
historical locks. The [selected snapshot](test_dependency_advisories_20261001.json)
retains all fifteen exact ranges and identifiers plus the complete raw-response
pin. No dismissal or historical-lock edit. The initial compact export failed
before writing because advisory objects lack `html_url`; the corrected export
constructs the standard GHSA URL from the actual identifier.

| Test Package | Previous | New |
| --- | ---: | ---: |
| Biopython | 1.86 | 1.87 |
| Pillow | 12.2.0 | 12.3.0 |
| setuptools | 81.0.0 | 83.0.0 |

New pins lie outside all fifteen retained ranges. The
[Biopython advisory](https://github.com/advisories/GHSA-x3vf-39hj-gxr4)
lists `<=1.86` with no explicit first patched version; range exclusion is not
comprehensive patch certification. The
[setuptools advisory](https://github.com/pypa/setuptools/security/advisories/GHSA-h35f-9h28-mq5c)
and [representative Pillow advisory](https://github.com/advisories/GHSA-9hw9-ch79-4vh6)
name 83.0.0 and 12.3.0 respectively. No exploitability or native/OS clearance.

Two new offline cases require current test pins to avoid all selected ranges
and prove the previous pins match all fifteen. NumPy stays 2.2.6. Application
requirements, scientific/default code, benchmark locks and scores are unchanged.

## Executed Validation

Fresh private Python 3.12 environment installation and `pip check` pass.
**511 focused cases pass in 71.77 seconds**, no failures/errors/skips, including
the two new cases, plot/count-reproduction modules, parser/reader checks and
packaging/archive/CPU-wheel orchestration. Twenty-two historical syntax warnings
remain. The [validation receipt](test_dependency_security_validation_20261001.json)
lists every executed module and pins JUnit, logs, inventory, artifacts and sources.
Counts overlap previous panels and are not additive.

The test pass initially uses ensurepip's pip 24.0. Patch pip to 26.2.1 privately
afterward; retain both inventories, separate installation receipt and a second
successful `pip check`. Do not relabel the 511-case pass as post-bootstrap-patch
execution. No shared interpreter or older private environment is modified.

These are checkout/retained-library checks, not fresh wheel installation,
whole-suite or full scientific benchmark equivalence. Post-push CI outcomes
and actual alert closure require separate observation. Historical vulnerable
locks remain evidence, not recommended fresh installations. Timing stays
deferred, without host polls/questions, DGX or unrelated job/service changes.
Uncertainty/source gaps, runtime/data rights and public release remain open;
the full publication goal is active and incomplete.
