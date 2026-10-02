# Standalone YGOB Numerical Reproduction

The frozen YGOB transfer analysis now has a standalone arithmetic path that
does not follow historical workstation files. The
[execution guide](../YGOB_ARITHMETIC_REPLAY.md) and
[compressed sufficient statistics](ygob_sufficient_counts_20261002.json.gz)
preserve the original four methods, 10,250 pillars, 12 point metrics and six
inferential endpoints. No inference, native admission, raw scoring, source
reference reconstruction, endpoint selection or parameter tuning is repeated.

## Findings and Scope

All 20,000 original paired multinomial replicates are regenerated using NumPy
2.2.6 / PCG64 seed 20260917. Integer half-count arithmetic reproduces every
effect and all 24 nominal/Bonferroni bound values within 1e-12. The primary
phylogenetic-OrthoHMM-minus-full-OrthoFinder F1 remains -0.084870 percentage
points, adjusted interval [-0.622528, 0.445222]. Reproduction does not establish
equivalence or superiority. Pillar exchangeability, cross-pillar FP allocation,
reference projection and homolog-family overlap remain explicit limitations.

The 77,438-byte gzip has SHA256
`62885bb85eb357d97e55d053134d262e1ebb1e2b088a37df96986f98f8620d66`.
Its 590,099 decoded bytes preserve 205,000 integer observation cells and shared
pillar sizes. Gene IDs, pillar labels, sequences, memberships and predictions
are absent; the original ordered pillar/size signature is retained only as a
hash. Projection checks the exact full-result and summary digests and all
shared metadata before export. A first `KeyError: full_results` attempt failed
before creating the output; the exporter now explicitly checks the exact
one-field summary/full-report schema difference. Original sources/results stay
unchanged. The first local numerical replay is retained at its earlier helper
identity, not relabeled as execution of the final stricter schema check.

## Relocated Archive

Build a local seven-member archive containing six selected payload files and
one externally pinned index: standalone script, counts, unchanged summary,
execution guide, pinned NumPy requirement and project license. Verify exact
archive/index identities, regular members and bounded inventory before fresh
extraction to `/tmp`; do not use `extractall` or trust an internal digest alone.
The actual fresh private Python 3.12.3 child runs with `-I -B` and no project
imports. It verifies every payload before and after reproduction. Three guard
canaries fail as intended; zero original-checkout, `/proc` or `/sys` opens occur
afterward. All six endpoints and 12 point metrics match. The child exits zero,
its copied report/logs remain local, and the temporary extraction is removed.

| Retained local artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `orthohmm-ygob-arithmetic-20261002-v1.tar.gz` | 89,006 | `794aa02ee7f91868eb2bc8492e92683d3c52afd9464d23eda66dbabdb24751f5` |
| `BUNDLE.json` inside archive | 1,589 | `242947e66b5e48d0c74533d0bd4b4cc7b5ac457ed8f2d5c6a591dfc16af6f1de` |
| Copied numerical replay report | 1,677 | `2c1a73572736549aaca5392a8df92913f48ddf4b0237f28f037a91f6320b2296` |

Archive location is `benchmarks/work/ygob_replay_bundle_20261002/`.
The [machine-readable receipt](ygob_arithmetic_replay_review_20261002.json)
pins the final script, all six selected source files, archive, restoration,
copied replay, original full result and tests. The one-off archive/audit helper
is retained locally with its own identity; it is not an installed public CLI.
The public standalone command needs just the script, counts and summary plus
the declared runtime; it does not need that one-off helper or archive.

Guards are Python-event checks, not OS containment. NumPy/runtime dependencies
load before guard installation and remain outside the archive. Same-host
restoration is not cross-platform or hermetic runtime proof. This local numerical
archive is not public deposition, raw-data rights clearance or the complete
executable publication release. Identifier removal alone does not clear an
arbitrary future bundle; raw YGOB acquisition and attribution limits remain.

## Regression and CI

The initial three-module panel has **50 passes in 1.59 seconds**, zero failures,
errors or skips. Its 32 new cases cover exact source/count identities, all point
metrics, integer half-counts, malformed scopes/counts/coverage/controls,
diagnostic exclusion, no-overwrite, bounded batching and wrong-digest failure
preservation. A small explicit resampling fixture independently checks ratios,
shared pairing and batch invariance. Interval-field tests use declared synthetic
draws rather than repeating the full 20,000-by-10,250 replay in every CI case.
The real full-size intervals are checked in the two retained actual executions.

After adding both manuscript availability links, the final four-module panel
has **62 passes in 1.75 seconds**, zero failures/errors/skips. It includes 34
new replay/prose cases, 18 existing YGOB scorer/report cases and ten historical
review-artifact cases. The two panels overlap and must not be added. Final
JUnit is 9,046 bytes, SHA256
`77c167cdc79bb4e05a1ffa135cbbb3b70b78ba82a9cac1aaee752da1d6ad94f6`.

At 13:23:44 UTC, preceding source `a485ba00` run 37010497782 has **all eight jobs
successful**. Download the terminal Python 3.10 job and its three-file artifact
once; actual checkout SHA, selection JSON and JUnit agree. Collection has 14,295
nodes/14,291 selected/four raw-source exclusions, zero collection errors/skips.
Unit execution has 14,167 passes/119 execution-time skips/four deselections;
all 25 profile cases pass. Its separate integration JUnit has five passes and
no skips/errors/failures. Sibling counts are not inferred. This confirms the
declared public profile, not the private raw/native gate or the newer YGOB code.

Main/extended availability text and reproduction/claim indexes link this
workflow. Earlier manuscript renders/archives remain their exact snapshots;
no new manuscript rendering or native run is performed. Controlled timing
stays deferred without host polls/questions, DGX or unrelated job/service
actions. Other-QfO uncertainty, complete runtime/rights/release/deposition and
final package reconciliation remain open. The original goal is incomplete.
