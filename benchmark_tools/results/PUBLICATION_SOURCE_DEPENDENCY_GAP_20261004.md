# Current Source Export Dependency Gap

A read-only inspection of committed workflow revision
`576d715a14e4535edd779a59607688832f1574ca` finds that the historical
`native-build` source profile selects 1,824 paths but omits 20 Python modules
under `benchmark_tools/results` imported by selected unit tests or their
transitive local dependencies. The [machine-readable observation](publication_source_dependency_gap_20261004.json)
records exact Git bytes/hashes and importing files. AST inspection traverses
1,834 committed Python modules without importing scientific code or running
native inference. This is static evidence, not complete runtime/data closure.
Independent readback verifies all 21 selector/dependency pins against both
the recorded Git revision and current bytes, every recorded direct import
edge, and actual selector exclusion for all 20 modules. No unit suite or
native analysis is rerun for this documentation-only observation.

Examples include the current table/figure/section generators, prepared
reporting bundler, terminal reviewer/capture and goal-reaffirmed launcher.
The historical source verifier compiles syntax without executing imports;
that bounded success cannot prove current exported tests collect or workflows
execute. Historical source archives and their verified scopes remain unchanged.
The separately executed reporting-only component remains valid within its
documented inventory; this finding does not invalidate it or live inference.

After the frozen timing panel is terminal and independently reviewed, extend
the final source export's explicit selection to include the operational
results-module sources and preserve them as exact committed Git blobs. Keep
historical profile/schema behavior separate. Test actual import/collection
from a relocated export as well as manifest verification; source syntax alone
is insufficient. The final handoff must also select the final manuscript,
resource table/figure/section and their real execution evidence. Do not label
the old source export as a complete current executable release.

The active timing recipe binds the top-level source bundler. It remains
unchanged during jobs 22422 and the remaining native identity; this finding
is publication packaging work, not a new timing gate, diagnostic rerun,
quiet-window prerequisite or reason to restart a native attempt. No unrelated
job or service is modified.
