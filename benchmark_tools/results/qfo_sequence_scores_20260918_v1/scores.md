# Corrected QfO Sequence-Search Controls

| Variant | Status | GO similarity | EC similarity | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | FAS | Secondary mean | Submitted pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_hits | admitted | 0.479327 | 0.874662 | 0.602031 | 0.628092 | 0.576306 | 0.715448 | 0.645978 | 11300151 |
| top100 | not_admitted | pending | pending | pending | pending | pending | pending | pending | pending |

- DIAMOND sequence-search controls with frozen OrthoHMM downstream grouping, not standalone competing tools.
- GO/EC similarity and FAS are not F1; the six-metric mean is project-defined and secondary.
- Native metric standard errors are not paired method-difference confidence intervals.
- Equal search cutoffs do not establish matched sensitivity or computational effort.
- The HMM baseline and paired uncertainty are not supplied by this table; no HMM advantage is established.
- Pending means no admitted score, not zero and not a scheduler-state assertion.
- Rechecks admitted hashes and native metric content, not the entire inference workflow.
