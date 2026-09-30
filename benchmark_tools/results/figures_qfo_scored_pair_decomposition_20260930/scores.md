# QfO Native Scored-Pair Decomposition

All contrasts are method minus full OrthoFinder 3.1.5. Values use serialized six-decimal pair scores.
Terms sum to the original rounded-mean difference; shared-pair fractions use each method's own eligible-pair count.
No new inference, scoring, confidence interval or causal attribution. GO/EC similarities are not F1.

| Metric | Method | Scored Pairs | Shared (% Method / Reference) | Difference (Points) | Shared Denominators | Method-Only | Negative Reference-Only |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GO | OrthoHMM sensitive | 145,619 | 53.81 / 47.91 | +0.272313 | +2.902999 | +20.757840 | -23.388525 |
| GO | OrthoHMM phylogenetic | 84,211 | 90.12 / 46.40 | +2.080130 | +21.571332 | +4.569645 | -24.060846 |
| GO | OrthoFinder sequence-only | 1,979,084 | 8.26 / 100.00 | -6.250487 | -43.074322 | +36.823836 | +0.000000 |
| GO | SonicParanoid | 197,930 | 59.17 / 71.61 | -1.516344 | -5.995204 | +16.911482 | -12.432622 |
| GO | ProteinOrtho | 89,609 | 87.33 / 47.85 | +1.678842 | +19.360431 | +5.812550 | -23.494138 |
| GO | FastOMA | 171,132 | 50.18 / 52.50 | -3.265446 | -1.127415 | +19.346573 | -21.484604 |
| GO | OrthoMCL | 175,181 | 64.38 / 68.95 | -0.579325 | -2.202827 | +15.380308 | -13.756805 |
| EC | OrthoHMM sensitive | 185,664 | 64.72 / 68.52 | -0.409252 | -3.680772 | +30.555622 | -27.284102 |
| EC | OrthoHMM phylogenetic | 117,460 | 86.59 / 58.00 | +2.952029 | +27.772815 | +12.451329 | -37.272115 |
| EC | OrthoFinder sequence-only | 1,660,940 | 10.56 / 100.00 | -18.248861 | -83.729418 | +65.480556 | +0.000000 |
| EC | SonicParanoid | 244,356 | 60.61 / 84.46 | -6.336924 | -22.795007 | +29.339215 | -12.881132 |
| EC | ProteinOrtho | 143,135 | 80.23 / 65.49 | +2.703790 | +14.350283 | +18.228315 | -29.874807 |
| EC | FastOMA | 178,598 | 62.57 / 63.72 | -3.896269 | -1.120597 | +29.009607 | -31.785278 |
| EC | OrthoMCL | 220,622 | 66.11 / 83.18 | -1.480422 | -16.404661 | +28.573735 | -13.649496 |
