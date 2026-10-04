# Actual Source And Review Components

Build actual components from committed workflow
`79a850e8779b808d0ac659ef563063909f3f989d`, rather than interpreting fixture
tests as final exports. The source includes result helpers under schema3 and
the `native-build` profile. Scientific source remains
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`; setup-only overlay remains
`6fd6df19daba83ec6467b917988f99e27a95be14`.

The [source index](final_publication_source_index_20261004.json) is 775,869
bytes, SHA-256
`e7f2f0e72eabbd45127f83ef75a9829c547bc611acc3034898349360672abf25`.
It binds 1,912 payload files: 43 scientific, 1,868 workflow and one build
overlay, totaling 11,828,926 bytes. All 1,891 Python files pass syntax checks.
The frozen parser/writer escape-sequence warnings remain recorded; no
scientific source is changed to silence them.

The [direct-review index](final_publication_review_index_20261004.json) is
30,609 bytes, SHA-256
`c6dcecc874d33e6ad06a404036c9a697356370607f07ac4abafdd21b64fff816`.
It binds 83 payload files totaling 6,577,344 bytes, with 57 direct targets and
59 HTML local-link occurrences. Explicit entrypoints select the 4 October
main text, HTML and twelve-page PDF. Its review/artifact revision is
`79a850e8`; ledger revision `904f4d4b` preserves the ledger actually rendered.
Historical default manuscript selections are not reinterpreted.

Both components are copied into fresh directories outside the repository.
Copied verifiers pass under isolated standard-library Python with
`PATH=/no-git`, using external manifest anchors and no original checkout imports.
All 20 previously omitted result helpers import from the copied workflow.
All 344 relevant cases in eight source/review/manuscript test modules collect;
project-module paths resolve beneath the copied workflow. These cases are
collected, not executed, and collection is not native inference reproduction.
The [execution receipt](final_publication_components_execution_20261004.json)
records actual build/verify/collection commands, outputs and canonical/copied
locations. Public indexes match canonical and copied bytes exactly.

Canonical components remain local under
`benchmarks/work/final_publication_components_20261004/{source,review}`.
No source/review payload archive is publicly uploaded here. Rights and private
runtime/data assets remain separate; a manifest is not redistribution clearance.
The direct-review component contains direct links only, not all transitive
evidence. Figures remain linked in its twelve-page main PDF.

## Remaining Work

Reconcile the complete figure appendix with the current twelve-page main text
and final resource figure; preserve old presentations and their exact inputs.
Build and validate the explicit-review handoff and executable/versioned archive
with the retained native runtime/data assets, final reporting component and
current scientific limitations. Record external release/deposition steps as
unexecuted unless actually performed. These component validations do not
establish whole-study restoration, isolated timing or publication readiness.
