# Canonical QfO Ordering Comparison

| Endpoint | Historical retained | Canonical | Canonical minus historical |
| --- | ---: | ---: | ---: |
| GO | 0.490349280 | 0.490349280 | +0.000000000 |
| EC | 0.965650500 | 0.965650500 | +0.000000000 |
| VGNC | 0.901690096 | 0.901690096 | +0.000000000 |
| SwissTrees | 0.833513221 | 0.833513221 | +0.000000000 |
| TreeFam-A | 0.614863704 | 0.614863704 | +0.000000000 |
| FAS | 0.762993312 | 0.762588129 | -0.000405183 |

Secondary six-metric means: 0.761510019 historical; 0.761442488 canonical. These are not F1 scores.

FAS uses unseeded native sampling; its difference is not an isolated ordering effect.
GO/EC are mean Schlicker scores; VGNC/SwissTrees/TreeFam-A use harmonic TPR/PPV.
Native axes, precision/recall and assessed-relation counts are retained in results.json.
NR_ORTHOLOGS is not necessarily total submitted-pair coverage.
No new paired confidence interval, independent generalization or default promotion.
