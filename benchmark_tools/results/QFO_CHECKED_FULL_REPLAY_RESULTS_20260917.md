# Checked Full QfO Replay

Audit-only job21480 completed0:0 in1:23 using frozen auditor
b38f3f5435dd9da7920884a10be741b9beabf8bb. The
[independent admission](qfo_checked_full_replay_recovered_verified_20260917.json)
has SHA256905c8accae2a0c98e73e028ec30e8a40d5a3d36d3223169afb6bedc01127758e
and status `checked_full_replay_recovered_verified`.

This recovers scientific output validation, not the failed wrapper's scheduler
status. Original inference job21333 remainsFAILED1:0 because its final guard
expected incorrect stage labels. The auditor accepted only the exact preserved
failure/report hashes and validated the actual frozen label contract. No
clustering, profile construction or other inference was rerun. The fresh input
check is explicitly retrospective because the wrapper failed before postflight.
See [failure and recovery details](QFO_REPLAY_LABEL_FAILURE_20260917.md).

## Validated Evidence

Allfour native clustering stages passed reconstructed graph, saved/native
fingerprint, execution environment, source, partition and provenance checks.
Allfour emitted partitions contain the full976504-gene universe without duplicate
or unknown memberships. There are379 checked provenance records.

| Output stage | Groups |
| --- | ---: |
| Multipass | 308190 |
| Multipass refined | 393231 |
| Strict profiles | 305656 |
| Strict profiles refined | 390980 |

The initial checked partition is byte-identical to the independently admitted
three-run initial-graph experiment:349898groups, SHA256
8c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd.
The replay built57883profiles and reports98085strict profile edges. These are
execution counts, not measures of accuracy or proof that those edges help.

## Historical Disagreement

The final historical partition has390817groups; this replay has390980.
Label/order-independent membership comparison finds4652historical-only groups
and4815replay-only groups. Both cover the same complete gene universe. Historical
scores therefore cannot be assigned to the recovered output, and the recovered
output must not silently replace the historical benchmark row.

Retain the checked replay as a separately identified baseline for subsequent
prespecified QfO analyses. This single full replay does not establish general
determinism, accuracy superiority, full historical reproducibility or matched
end-to-end efficiency. Accuracy, uncertainty and HMM/phylogeny contribution
analyses remain due; no accuracy score was used to select this partition.
