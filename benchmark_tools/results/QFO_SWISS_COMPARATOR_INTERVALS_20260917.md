# SwissTrees Comparator Paired Intervals

18 families; 100000 shared draws; seed20260920. Candidate minus reference in raw0-to1 units.
Adjustment covers24endpoints. Historical development-exposed comparison, not independent confirmation.

| Candidate | Reference | Metric | Difference | Nominal95% CI | Adjusted CI | Wins/ties/losses |
| --- | --- | --- | ---: | --- | --- | --- |
| orthohmm_high_sensitivity | orthofinder_3_1_5_full | F1 | -0.185310 | [-0.248508, -0.126842] | [-0.285081, -0.091574] | 1/0/17 |
| orthohmm_high_sensitivity | orthofinder_3_1_5_full | PPV | -0.299852 | [-0.391424, -0.206409] | [-0.438896, -0.156669] | 1/0/17 |
| orthohmm_high_sensitivity | orthofinder_3_1_5_full | TPR | -0.082121 | [-0.135504, -0.026894] | [-0.162660, 0.006846] | 2/5/11 |
| orthohmm_phylogeny_satellite_v2 | orthofinder_3_1_5_full | F1 | -0.067184 | [-0.104124, -0.035760] | [-0.127608, -0.021232] | 0/5/13 |
| orthohmm_phylogeny_satellite_v2 | orthofinder_3_1_5_full | PPV | +0.004959 | [-0.019786, 0.041027] | [-0.029789, 0.068449] | 3/5/10 |
| orthohmm_phylogeny_satellite_v2 | orthofinder_3_1_5_full | TPR | -0.108669 | [-0.155406, -0.062996] | [-0.182107, -0.040847] | 0/5/13 |
| orthofinder_3_1_5_sequence_only | orthofinder_3_1_5_full | F1 | -0.159653 | [-0.260259, -0.069788] | [-0.318956, -0.022866] | 2/2/14 |
| orthofinder_3_1_5_sequence_only | orthofinder_3_1_5_full | PPV | -0.345148 | [-0.464114, -0.222921] | [-0.524955, -0.153102] | 2/2/14 |
| orthofinder_3_1_5_sequence_only | orthofinder_3_1_5_full | TPR | +0.051934 | [0.008793, 0.106380] | [0.001029, 0.143955] | 6/12/0 |
| sonicparanoid_2_0_9 | orthofinder_3_1_5_full | F1 | -0.082054 | [-0.125170, -0.045126] | [-0.151991, -0.029443] | 0/1/17 |
| sonicparanoid_2_0_9 | orthofinder_3_1_5_full | PPV | -0.077697 | [-0.140216, -0.015551] | [-0.177787, 0.019127] | 1/1/16 |
| sonicparanoid_2_0_9 | orthofinder_3_1_5_full | TPR | -0.083733 | [-0.130209, -0.042852] | [-0.158180, -0.025828] | 0/3/15 |
| proteinortho_6_3_6 | orthofinder_3_1_5_full | F1 | -0.161325 | [-0.226319, -0.104004] | [-0.264228, -0.075970] | 0/2/16 |
| proteinortho_6_3_6 | orthofinder_3_1_5_full | PPV | +0.006035 | [-0.013718, 0.037455] | [-0.018873, 0.063292] | 2/3/13 |
| proteinortho_6_3_6 | orthofinder_3_1_5_full | TPR | -0.236700 | [-0.315519, -0.160654] | [-0.362345, -0.119480] | 0/2/16 |
| fastoma_0_3_5 | orthofinder_3_1_5_full | F1 | -0.096433 | [-0.140754, -0.060756] | [-0.171393, -0.043898] | 1/0/17 |
| fastoma_0_3_5 | orthofinder_3_1_5_full | PPV | -0.037437 | [-0.101296, 0.019951] | [-0.142082, 0.052141] | 4/0/14 |
| fastoma_0_3_5 | orthofinder_3_1_5_full | TPR | -0.129950 | [-0.195095, -0.072178] | [-0.237141, -0.045494] | 1/3/14 |
| orthomcl_1_4 | orthofinder_3_1_5_full | F1 | -0.100802 | [-0.148940, -0.058703] | [-0.177237, -0.038303] | 1/2/15 |
| orthomcl_1_4 | orthofinder_3_1_5_full | PPV | -0.043938 | [-0.100215, 0.009761] | [-0.135083, 0.039583] | 2/2/14 |
| orthomcl_1_4 | orthofinder_3_1_5_full | TPR | -0.133032 | [-0.202621, -0.069019] | [-0.242410, -0.038786] | 1/6/11 |
| orthohmm_phylogeny_satellite_v2 | orthohmm_high_sensitivity | F1 | +0.118125 | [0.039714, 0.197221] | [-0.007676, 0.241518] | 15/1/2 |
| orthohmm_phylogeny_satellite_v2 | orthohmm_high_sensitivity | PPV | +0.304810 | [0.205414, 0.403840] | [0.150217, 0.455082] | 16/1/1 |
| orthohmm_phylogeny_satellite_v2 | orthohmm_high_sensitivity | TPR | -0.026547 | [-0.085657, 0.017477] | [-0.128297, 0.031889] | 5/8/5 |

- Retrospective, development-exposed; not independent or model-selection-adjusted confirmation.
- Only 18 curated families; shared history and merged predictions can violate exchangeability.
- Approximate conditional family-bootstrap sensitivity analysis; adjustment covers these 24 endpoints only.
- Family wins/ties/losses are descriptive; inclusion of zero does not establish equivalence.
- Sequence-only OrthoFinder is an MCL-checkpoint diagnostic; FastOMA used a supplied OrthoFinder tree.
- Phylogenetic versus sensitive OrthoHMM is not a pure reconciliation ablation.
- No interval transfer to recovered ablations, other QfO challenges or the secondary six-metric mean.
