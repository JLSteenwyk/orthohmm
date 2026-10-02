# Retained Native FAS Test Reference

## Exact Source, Explicit Selection

Two FAS population-audit tests previously read the external
`qfo_benchmark/benchmark-webservice/fas_benchmark.py`. That checkout is absent
from macOS CI; the actual source-00dd Python 3.10 fast log reports both failures.
Do not replace the native function with an invented stub or skip the comparison.

Add the complete unchanged source and upstream license under
`tests/samples/native_qfo_fas`, with a provenance inventory and scope notice.
Both match the retained upstream Git commit
`c0854a96c1a0fd7f2a891d971af0863002fabc90` and independently fetched public
[source](https://raw.githubusercontent.com/qfo/benchmark-webservice/c0854a96c1a0fd7f2a891d971af0863002fabc90/fas_benchmark.py)
and [license](https://raw.githubusercontent.com/qfo/benchmark-webservice/c0854a96c1a0fd7f2a891d971af0863002fabc90/LICENSE)
blobs. The source also exactly matches the original canonical source archive
member, whose whole-archive SHA256 was checked before comparison. No upstream
source modifications, data, annotations or Darwin binary are included. The
full upstream MPL 2.0 license and separate Darwin notice remain attached;
these files do not acquire the OrthoHMM project license.
Their existing trailing whitespace and license final blank line are preserved
to retain exact upstream hashes; staged whitespace checks pass for owned edits
when only those two immutable upstream files are excluded.

Tests explicitly select this fixed reference, not an automatic fallback.
Independent hardcoded size/hash pins check the source and license before AST
use. The native loader function alone runs on synthetic lookup entries; native
overwrite/invalid-value behavior and the exact query literal are still compared
with the actual population auditor. Module-level external imports and the
full FAS benchmark are not executed. Six altered/missing/symlink cases fail
closed; a seventh new case checks provenance and scope.

## Actual Validation

**65 local cases pass in 0.84s**, zero failures/errors/skips: 33 FAS cases plus
the preceding 32 provider cases. In a fresh copied tree outside the checkout,
**all 33 FAS cases pass in 0.68s**. The two panels overlap, not 98 independent
cases. Shared installed dependencies remain in use.

The child uses isolated Python flags, a copy-first import path, disabled pytest
plugin autoload and a Python audit hook forbidding subprocesses and original
checkout opens/chdir. The canary is rejected, subsequent original-path events
are zero, all 25 loaded project modules originate in the copy, and the full
native scorer module is not imported. The copy contains 877 files: 872 committed
Python files and five explicitly pinned owned overlays, totaling 6,078,966 bytes, and no
external QfO checkout. Temporary staging is removed after successful execution.
[The compact copied receipt](fas_native_source_copied_receipt_20261002.json)
retains the exact guard, invocation, observations, JUnit identity and full local
staging-receipt hash. This is same-host Python-level guarding, not OS containment
or complete independent/native scoring reproduction.

[Machine-readable evidence pins](fas_native_source_reference_20261002.json)
also bind local tests, source/license/provenance/notice bytes and the inspected
prior CI log. At 06:46:05 UTC source-657 automatic run 36974478584 has five live
test jobs; wheel/docs pass. Source-00dd is terminal with all five test jobs
failed, wheel/docs successful. Only its already-retained Python 3.10 fast log
is used for the source-availability diagnosis. No new-FAS macOS confirmation
or sibling/full-suite success is claimed; no CI handle is restarted.

Later inspect actual source-657 Python 3.10 fast CI: all 138 preceding
provider/source/pair-count/checker/bundler cases pass, including all 32 changed
provider cases. Overall: 13,857 passes, 29 failures, 110 skips, 30 warnings,
505.47s. The two old external FAS-source failures remain. The old provider
failure disappears, but owned-command timeout has a new EPERM failure; the
unchanged overall failure count is not evidence that the provider fix failed
or that timeout handling is fully resolved. New FAS tests are not in this
source. At 06:50:35 UTC full/3.10/3.12/3.13 test jobs fail, 3.11 remains live,
wheel/docs succeed. Only the 3.10 log is inspected.

## Boundaries

Production population-audit logic, native source-admission rules, frozen method,
scientific scores, source paths and historical receipts stay unchanged. Earlier
59.96-million-lookup population analysis and its negative findings are not rerun
or relabeled. The snapshot is not a complete QfO installation, upstream data
rights clearance, a remedy for the unseeded FAS sample, an uncertainty estimate,
controlled timing, a complete executable release or public deposition.
