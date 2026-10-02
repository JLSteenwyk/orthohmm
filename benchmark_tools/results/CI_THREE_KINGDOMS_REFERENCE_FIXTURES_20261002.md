# Three Kingdoms Assessment Fixture Boundaries

The matched SonicParanoid assessment test must exercise real conversion,
scoring and independent pair-count comparison without silently assuming that
untracked biological data ships with a fresh checkout. This test-only change
adds an explicitly synthetic reference branch and retains a separately labeled,
hash-checked branch for the original BUSCO reference. It changes neither the
production reference nor any scientific score or admission rule.

## Observed Failure

The retained source-42e macOS Python 3.11 fast log shows four integrated fixture
errors copying `three_kingdoms/busco/reference_orthogroups.txt`, which is absent
from the checkout. Reuse that pinned log; do not infer new CI test outcomes.
At 03:49:06 UTC on 2 October the existing source-d94 CI has five live test jobs,
wheel success and queued docs. Its Python 3.11 fast job API reports successful
native-utility setup, but no new log, utility versions or test outcomes have
been inspected. Setup-step success does not establish regression success.

Later, source-d94 Python 3.10 fast job `110695664741` completes failure: 13,701
passed, 55 failed, four errors, 105 skipped and 30 warnings in 463.80 seconds.
Its actual log records GNU Time 1.10 and Pandoc 3.11. Within the prior 100-case
utility scope, 93 pass and seven fail without skips/errors; native GNU Time,
parser, search and simulation cases pass. The seven remaining QfO runner cases
return exit 127 when their fixtures invoke `/bin/true` or `/bin/false`. Their
native-log contents are not inspected, so the exact child launch cause remains
unverified. The four missing-reference fixture errors remain in that old source.
At 03:56:49 UTC all four fast jobs fail, the full job remains live and wheel/docs
succeed. Only the Python 3.10 fast log is newly inspected; no sibling cause or
full-suite success inference follows.

## Two Explicit Reference Branches

The portable branch generates distinct `synthetic_family_*` identifiers with
255 groups, 2,035 genes and 7,352 within-reference gene pairs. Its arbitrary
size histogram is 238 groups of eight genes, ten of seven, four of six, two of
four and one of 29. A separate arithmetic test verifies these dimensions and
uniqueness. The dimensions intentionally match the unchanged production
universe guard; they do not reconstruct BUSCO families, taxonomic membership,
identifiers or the retained family-size distribution. No biological reference
is inferred or substituted in a real benchmark.

The retained-data branch checks the actual reference bytes/SHA-256 against the
original committed SonicParanoid plan before copying it. Missing data is an
explicit skip for that branch only. Available but changed data fails; it is
never replaced with synthetic truth. The retained local file is unchanged at
30,699 bytes, SHA-256
`a5f3447056ecfa305442caff0d898d524eed350f13587d3247d3e1ba4c19757d`.

Both branches exercise the same actual conversion and scoring subprocesses,
native-table/input consistency checks and independent integer pair counter.
Inference, scheduler identity and runtime gates remain synthetic fixtures;
there is no real SonicParanoid run or new scientific accuracy admission.
Perfect predictions yield exactly 7,352 TP, zero FP/FN and F1 1. The separated
member case now checks exact lost-pair arithmetic and F1, not just nonzero FN.
Corrupt conversion/scoring output remains rejected in both branches.

Two new guard cases in each branch reject changed reference bytes before an
assessment directory is created, and reject a changed universe even when the
temporary fixture's plan/hash binding is explicitly updated. The latter uses
the unchanged production dimension guard after real scoring; failed evidence
is retained with `accuracy_admitted: false`. No production repinning authority
or relaxation is introduced.

## Portable QfO Smoke Commands

Replace only the three QfO runner fixtures' tiny true/false commands with
`sys.executable -c` and explicit `SystemExit(0)` or `SystemExit(1)`. This preserves
real child execution and success/failure assertions without requiring a specific
location for external true/false utilities. The production commands and GNU Time
binding are unchanged. Existing allocation, provenance, pending-admission,
runtime-drift and no-restart checks still run; no exit-127 acceptance is added.
New remote confirmation is required rather than inferred from local success.

## Executed Validation

**58 focused Linux cases pass in 1.87 seconds**, zero failures/errors/skips,
including nine additional parameterized cases. Module counts: assessment 34,
independent pair-count audit 15, normalizer five and scorer four. Both retained
and synthetic integration branches execute locally.

In a separate process the retained source path is deliberately reported absent
before test-module import, without deleting or editing it. The same scope gives
**52 passes and six explicit retained-data skips in 1.36 seconds**, zero
failures/errors. All six synthetic integration cases still execute. This is a
local capability simulation, not actual macOS; the panels are not additive.

The initial focused run had two assertion failures because new tests expected
the wrong exception string, while production correctly rejected changed bytes.
Only the test expectation was corrected. A later exact-pair strengthening
supersedes the earlier passing run; all earlier reports remain pinned rather
than discarded. All reports remain recorded in the
[machine-readable receipt](ci_three_kingdoms_reference_fixtures_20261002.json).
All six production assessment/conversion/scoring/validation modules and the
historical plan are byte-equal to the base commit. The original reference and
every retained scientific endpoint stay unchanged.

The combined utility/parser and Three Kingdoms scope, after the three QfO
fixture-command changes, gives **158 passes in 5.12 seconds**, zero failures,
errors or skips. A separate combined process with retained data absent gives
**152 passes and six retained-data skips in 4.69 seconds**, zero failures/errors.
These combined reports include the 58-case scope; do not add their observations
to the earlier panels. Three unchanged QfO production runners and three changed
test source pins are included alongside the two new combined JUnit reports.

## Remaining Requirements

Commit/push and inspect actual new CI without manually restarting existing jobs.
New remote fixture outcomes, remaining workstation-path/utility/Linux-specific
failures and broader executable-release work are open. These fixtures do not
redistribute BUSCO data, clear data rights, establish biological accuracy or
admit controlled runtime/memory measurements. Real reference restoration still
requires the correct external artifact and its provenance.

No shared package upgrade, native scientific rerun, local workload poll/question,
DGX work or unrelated job/service action occurs. Controlled timing remains
deferred. Remaining QfO uncertainty, source/rights review, comparable resources
and the final versioned release/archive keep the full publication goal active.
