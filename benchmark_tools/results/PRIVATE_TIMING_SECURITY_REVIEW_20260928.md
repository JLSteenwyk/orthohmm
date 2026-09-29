# Private Timing Security Review

The [retained repository snapshot](dependency_alerts_private_timing_20260928.json)
contains 26 open alerts, compared with 12 in the preceding release review.
The 14 added alert numbers (71-84) all refer to the first private timing lock:
13 Pillow advisories and one additional manifest occurrence of the existing
setuptools advisory. No alert was dismissed and no historical lock changed.
The absence of alerts naming the v2 lock is not evidence that its pins are safe.

The [live candidate audit](private_timing_advisory_ranges_20260928.json)
checks all 30 installed distribution versions against the retained aligned
candidate selection before comparing advisory ranges. That inventory matches.
The evaluator is shared with the prior release audits; comparisons distinguish
an absent package from an installed version outside a reported range.

| Installed package | Version | Matching unique advisories | Matching manifest alerts |
|---|---|---:|---:|
| Pillow | 12.2.0 | 13 | 13 |
| setuptools | 81.0.0 | 1 | 2 |
| Total | | 14 | 15 |

These are version-range matches, not confirmed exploitable inference paths.
Neither the current candidate nor the previous zero-match reader/recovery
checks have comprehensive security clearance. The new Pillow advisories were
not part of the earlier 12-alert comparison.

## Next Deployment Revision

Prepare a separate, hash-pinned private candidate using Pillow 12.3.0 and
setuptools 83.0.0, retaining all other selected versions and the frozen
scientific source. These are the first patched versions listed in the snapshot.
The [Pillow release notes](https://pillow.readthedocs.io/en/stable/releasenotes/12.3.0.html)
describe security fixes as well as other changes; do not assume output or
runtime equivalence. The [setuptools advisory](https://github.com/pypa/setuptools/security/advisories/GHSA-h35f-9h28-mq5c)
concerns source-distribution exclusions on macOS, not an established exploit
in this Linux inference pipeline.

Preserve the existing candidate and its installation receipts. Validate the
new full inventory, dependency closure, imported payload changes and native
fixture predictions before updating the timing execution amendment. Recheck
advisory ranges without presenting that bounded check as full security review.
No upgrade or scientific run was performed in this audit. YAML build provenance,
collector integration, final resource accounting and quiet-window verification
remain separate timing requirements.

Thirteen focused tests pass, including duplicate advisory accounting, distinct
absent/patched statuses, inventory drift and malformed snapshot scope rejection.
Slurm returned an empty queue during this turn; this is not a whole-host quiet
assessment and no production timing was launched. No DGX access occurred.
