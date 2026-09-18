# Matched-Input SonicParanoid Result

Inference21795 completed0:0 in1:17:55 onbizon with32CPUs/192GiB.
Assessment21796 completed0:0 in1:24. The independent assessment status is
`matched_three_kingdoms_sonic_score_verified`, accuracy_admitted=true.

| Metric | Value |
|---|---:|
| Pair F1 | 0.9912758996728462 |
| Precision | 0.9934426229508196 |
| Recall | 0.9891186071817193 |
| True-positive pairs | 7272 |
| False-positive pairs | 48 |
| False-negative pairs | 80 |
| Reference gene coverage | 2031/2035 |
| Predicted groups | 19870 |
| BUSCO reference families | 255 |

The run uses the frozen matched12-proteome input, resolving the historical
Sonic input mismatch for this contemporary comparison. It does not overwrite
the historical result or isolate a causal effect of the seven affected Danio
sequences. All retained checked records, normalized output, conversion log,
scorer output and source were rehashed; F1 was independently recomputed as
2TP/(2TP+FP+FN).36focused runner/assessment tests pass.

[Assessment report](three_kingdoms_sonic_matched_assessment_21796.json),
SHA-256`aeddd422a5c733ac40bef0c4953a6bda12ab2a751baeeae9e9c28cac472e9653`.
Frozen inference executor7c3d0784d1a22da0a0e37db6860977c6649fd2ca;
assessment executor729af62e5d6cf9e449eafe16190b2a01062201ec.

This supplementary endpoint measures BUSCO reference-gene co-membership.
It does not penalize predictions outside the reference, establish genome-wide
accuracy, or test whether phylogenetic pairs are better than group cliques.
Shared-host elapsed time is descriptive, not controlled scaling evidence.
Updated aggregate figures/tables must distinguish this matched-input result
from the retained historical Sonic value rather than silently mixing them.
