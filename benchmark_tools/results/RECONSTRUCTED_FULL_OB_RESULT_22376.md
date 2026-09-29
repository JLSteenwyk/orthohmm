# Full Reconstructed-Base OrthoBench Reproduction

Job 22376 and its batch step completed with exit 0:0 on `bizon`, using
32 CPUs and 128 GiB. Scheduler elapsed time was 02:53:26 on the shared host;
this is execution provenance, not a controlled timing measurement.

The [presubmission protocol](RECONSTRUCTED_FULL_OB_PROTOCOL_20260929.md)
changed only the controller/base/installer Python and output path relative
to job 22337. No scientific parameter changed, checkpoint was reused, or
failed attempt retried. All eight workflow stages completed with their
prespecified commands and environments.

## Independent Verification

[Machine-readable result](reconstructed_full_ob_result_22376.json)
(39,949 bytes; SHA256
`81c4de02279aa08500c9275a611cf180a48f64c3891234848a1c6a0f81da0fc0`)
records the separately executed admission. It checked the pinned plan and
inputs, both installed package inventories and wheel payloads, and four
independent structure/sequence/event/hierarchy readers. The successful
full-data execution of this checker follows its focused unit tests.

| Endpoint | Reconstructed-base result | Comparison with job 22337 |
| --- | ---: | --- |
| Input genes | 251,378 | Same universe |
| Root groups | 59,770 | All identical, ignoring labels |
| Reference families | 70 | Complete score objects identical |
| Precision | 81.77045380181866% | No change |
| Recall | 67.75533630827974% | No change |
| F1 | 74.10607351873405% | No change |
| Native ortholog pairs | 966,439 | Pair file byte-identical |
| Reconciled candidate families | 8,681 | Same count |
| Bypassed candidate families | 45,764 | Same count |

Pair-confidence, reconciliation-node and hierarchical-group TSVs are also
byte-identical to the admitted baseline. Compared output checksums are bound
to independent scientific readbacks on both sides. There are no changed root
groups or reference-family score records. Both historical and new results
are retained; no score row is silently replaced.

## Scope

This establishes full-data reproduction on the same host using the separately
reconstructed base Python and fresh inference/reader environments. It closes
the earlier fixture-only evidence gap for that execution path. Existing
native assets, readers and inputs were supplied locally; this is not complete
archive restoration, OS isolation, cross-host validation, a new biological
test, or redistribution/security clearance. Wheel payload exclusions remain
explicit. The 27-run timing panel, remaining scientific uncertainty and
publication/release requirements are not completed by this result.
