# Conditional Native FAS Design Ranges

Observed Z is preserved separately from the simultaneous conditional 95% range
for expected repeated native ratio theta. Unknown G and mu_G preclude an exact
theta point value. These are not biological error bars on Z or unconditional
historical certificates. Table decimals are rounded; JSON retains full values.

| Method | Observed Z | Conditional Theta Range | G Bounds |
| --- | ---: | --- | --- |
| orthohmm_high_sensitivity | 0.776585042 | [0.771332967, 0.780286193] | [2125583, 2129667] |
| orthohmm_phylogeny_satellite_v2 | 0.762993312 | [0.761929649, 0.763906702] | [304262, 305332] |
| orthofinder_3_1_5_full | 0.691422133 | [0.688429518, 0.694915139] | [2384608, 2392806] |
| orthofinder_3_1_5_sequence_only | 0.561753084 | [0.546499945, 0.577724179] | [133118353, 133664208] |
| sonicparanoid_2_0_9 | 0.736680187 | [0.732115268, 0.738512044] | [2568174, 2572794] |
| proteinortho_6_3_6 | 0.813594864 | [0.813565174, 0.813689528] | [15385, 15394] |
| fastoma_0_3_5 | 0.654366336 | [0.643583534, 0.662539506] | [7331620, 7366500] |
| orthomcl_1_4 | 0.733752216 | [0.731706501, 0.735008398] | [1063057, 1090063] |

All 28 contrasts use left-minus-right in the frozen method order.

| Left | Right | Observed Z Difference | Conditional Expected Difference |
| --- | --- | ---: | --- |
| orthohmm_high_sensitivity | orthohmm_phylogeny_satellite_v2 | 0.013591730 | [0.007426264, 0.018356544] |
| orthohmm_high_sensitivity | orthofinder_3_1_5_full | 0.085162909 | [0.076417827, 0.091856675] |
| orthohmm_high_sensitivity | orthofinder_3_1_5_sequence_only | 0.214831958 | [0.193608787, 0.233786248] |
| orthohmm_high_sensitivity | sonicparanoid_2_0_9 | 0.039904855 | [0.032820923, 0.048170926] |
| orthohmm_high_sensitivity | proteinortho_6_3_6 | -0.037009822 | [-0.042356561, -0.033278981] |
| orthohmm_high_sensitivity | fastoma_0_3_5 | 0.122218705 | [0.108793461, 0.136702659] |
| orthohmm_high_sensitivity | orthomcl_1_4 | 0.042832826 | [0.036324569, 0.048579692] |
| orthohmm_phylogeny_satellite_v2 | orthofinder_3_1_5_full | 0.071571179 | [0.067014510, 0.075477184] |
| orthohmm_phylogeny_satellite_v2 | orthofinder_3_1_5_sequence_only | 0.201240228 | [0.184205470, 0.217406757] |
| orthohmm_phylogeny_satellite_v2 | sonicparanoid_2_0_9 | 0.026313125 | [0.023417605, 0.031791434] |
| orthohmm_phylogeny_satellite_v2 | proteinortho_6_3_6 | -0.050601552 | [-0.051759879, -0.049658472] |
| orthohmm_phylogeny_satellite_v2 | fastoma_0_3_5 | 0.108626976 | [0.099390143, 0.120323168] |
| orthohmm_phylogeny_satellite_v2 | orthomcl_1_4 | 0.029241097 | [0.026921251, 0.032200201] |
| orthofinder_3_1_5_full | orthofinder_3_1_5_sequence_only | 0.129669049 | [0.110705339, 0.148415194] |
| orthofinder_3_1_5_full | sonicparanoid_2_0_9 | -0.045258054 | [-0.050082526, -0.037200128] |
| orthofinder_3_1_5_full | proteinortho_6_3_6 | -0.122172731 | [-0.125260009, -0.118650035] |
| orthofinder_3_1_5_full | fastoma_0_3_5 | 0.037055797 | [0.025890012, 0.051331605] |
| orthofinder_3_1_5_full | orthomcl_1_4 | -0.042330082 | [-0.046578880, -0.036791362] |
| orthofinder_3_1_5_sequence_only | sonicparanoid_2_0_9 | -0.174927103 | [-0.192012099, -0.154391088] |
| orthofinder_3_1_5_sequence_only | proteinortho_6_3_6 | -0.251841780 | [-0.267189583, -0.235840994] |
| orthofinder_3_1_5_sequence_only | fastoma_0_3_5 | -0.092613253 | [-0.116039561, -0.065859355] |
| orthofinder_3_1_5_sequence_only | orthomcl_1_4 | -0.171999132 | [-0.188508453, -0.153982322] |
| sonicparanoid_2_0_9 | proteinortho_6_3_6 | -0.076914677 | [-0.081574260, -0.075053130] |
| sonicparanoid_2_0_9 | fastoma_0_3_5 | 0.082313851 | [0.069575762, 0.094928510] |
| sonicparanoid_2_0_9 | orthomcl_1_4 | 0.002927972 | [-0.002893130, 0.006805543] |
| proteinortho_6_3_6 | fastoma_0_3_5 | 0.159228528 | [0.151025668, 0.170105993] |
| proteinortho_6_3_6 | orthomcl_1_4 | 0.079842649 | [0.078556776, 0.081983026] |
| fastoma_0_3_5 | orthomcl_1_4 | -0.079385879 | [-0.091424864, -0.069166995] |

## Scope

- Conditional uniform fixed-population sampling and fixed pair outcomes are assumptions, not historical certification.
- Observed fixed means are separate descriptive scores, not point values or centers of the expected-ratio ranges.
- Known precomputed population means are conditioned on the retained audit; historical parser/database hash gaps persist.
- Finite context checks do not establish every historical annotation, file race, batch failure or worker choice.
- No biological pair/family or cross-method independence, missing-at-random claim, new score, sample or rescore.
- No interval for another endpoint, overall tool superiority, independent biological confirmation or publication readiness.
