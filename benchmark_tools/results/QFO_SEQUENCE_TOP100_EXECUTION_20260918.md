# Top-100 Sequence Graph Execution

Graph21814 completed0:0 in35:24. Its execution report status is
`sequence_checked_graph_complete_pending_admission`, not accuracy admission.
The984137-gene partitions contain217371multipass and414510multipass_refined
groups. Both outputs are23,875,927bytes. The refined output SHA-256 is
`db1fcf5d9c9dd364afeebf1ad7455c5951cc40613aad7bf25337c607b9b4a86e`.

[Execution report](qfo_sequence_top100_execution_21814.json), SHA-256
`ff401cf6936c3a4eafdb953f421b9008a931c127b0c022a7112439f680d32ba9`.
Plan SHA-256 remains
`0e20d71936b92d7d2f02f3f35e7b7f380553803d9e482fde4d4b6c549b6e96d0`.

Verified the clean frozen independent validator
f1e21b09c28f270dc3ef2243bdcad86f212b58a0 and submitted21827 with the
exact report hash. Scheduler confirms RUNNING2CPUs/192GiB. Pair conversion
must wait for successful validation. The diagnostic top100arm does not
replace the all-hit arm; no parameter or endpoint was changed after results.
Elapsed time is shared-host execution evidence, not matched scaling.

## Independent Admission And Downstream Jobs

Admission21827 subsequently completed0:0 in10:49 with status
`corrected_sequence_graph_admitted`. [Retained report](qfo_sequence_graph_admission_top100_21827.json)
SHA-256:`5477d4ac8288d0362619b3ce7613b8f12de722b0ba5cf767a3c7ab60b2bb9905`.
No accuracy endpoint was evaluated at this gate.

After verifying the clean frozen converter4f0c30e5cdf287a35c9600886aec0a41bcc0b720,
submitted pair conversion21828 with the exact admission identity; confirmed
RUNNING2CPUs/64GiB. Queued native scoring21829afterany:21828 using frozen
425f0a7f5d5dc9e1438ab0a1766c45b13596dc14 and independent admission21830
afterany:21829 usinge798c159b9f65727d5b1055b98a81cba3a8f374d. These match
the all-hit arm's gated workflow. Failed upstream stages cannot pass the
downstream provenance/completion checks. No score is admitted yet.
