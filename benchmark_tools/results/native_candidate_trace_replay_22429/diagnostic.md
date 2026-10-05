# Accepted Candidate Trace Replay

Job 22429, p0_c1_r0; read-only retained-data diagnostic.

Both complete partitions reconstruct from 63,245 original seeds and accepted unions.
Each has 54,745 groups covering 251,378 genes; 58 genes occur in changed groups.
Common accepted merges: 8,482 / 8,500.

| Round | Anchor gene | Old / native attachments | Both at cap | Changed-selection support spread |
| ---: | --- | --- | --- | ---: |
| 0 | WBGene00001875.1 | 4 / 4 | True | 0 |
| 0 | WBGene00001915.1 | 4 / 4 | True | 7.1054273576010019e-15 |
| 0 | WBGene00001877.1 | 4 / 4 | True | 0 |

Complete changed memberships, selections, numeric differences and direct pins are in diagnostic.json.

- Replay of accepted unions reconstructs whole partitions; candidate eligibility, rejected alternatives and score computation are not replayed.
- The retained original seed is sufficient for both traces. This does not independently establish every fresh seed's numeric order or the provenance of score-bit differences.
- Attachment caps and near ties describe observed decisions, not a demonstrated causal numerical mechanism or a reason to change the frozen method.
- Reference exclusions are inherited from the pinned separate score, not new independent accuracy evidence.
- Shared-host timing distortion remains unknown and potentially tool-dependent; this diagnostic does not attribute membership variation to contention.
