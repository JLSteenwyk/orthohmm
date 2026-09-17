# Recovered QfO Stage Results

Descriptive results; paired uncertainty remains outstanding.

| Stage | GO | EC | VGNC | SwissTrees | TreeFam-A | FAS | Secondary mean | Retained pairs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| multipass | 0.444803 | 0.884110 | 0.276207 | 0.682896 | 0.679114 | 0.773707 | 0.623473 | 32012674 |
| multipass_refined | 0.472537 | 0.931175 | 0.668209 | 0.677567 | 0.576755 | 0.775253 | 0.683583 | 8710340 |
| strict_profiles | 0.444006 | 0.884420 | 0.275058 | 0.679189 | 0.678899 | 0.771575 | 0.622191 | 31682321 |
| strict_profiles_refined | 0.472350 | 0.931347 | 0.667695 | 0.673851 | 0.576358 | 0.774633 | 0.682706 | 8710722 |

## Prespecified Differences

Candidate minus reference; raw metric units, not percentage changes.

| Contrast | GO | EC | VGNC | SwissTrees | TreeFam-A | FAS | Secondary mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| strict_profiles - multipass | -0.000797 | +0.000309 | -0.001149 | -0.003706 | -0.000216 | -0.002132 | -0.001282 |
| strict_profiles_refined - multipass_refined | -0.000186 | +0.000173 | -0.000514 | -0.003716 | -0.000398 | -0.000619 | -0.000877 |
| multipass_refined - multipass | +0.027734 | +0.047064 | +0.392002 | -0.005329 | -0.102359 | +0.001546 | +0.060110 |
| strict_profiles_refined - strict_profiles | +0.028344 | +0.046928 | +0.392637 | -0.005338 | -0.102541 | +0.003058 | +0.060515 |

## Native Coordinates

Native challenge-specific axes and counts are retained, not pooled across challenges.

| Stage | Challenge | X axis | X | Y axis | Y |
| --- | --- | --- | ---: | --- | ---: |
| multipass | GO | NR_ORTHOLOGS | 353317 | avg Schlicker | 0.44480321 |
| multipass | EC | NR_ORTHOLOGS | 313327 | avg Schlicker | 0.88411008 |
| multipass | VGNC | TPR | 0.9689980780479652 | PPV | 0.16105779246933985 |
| multipass | SwissTrees | TPR | 0.72172476 | PPV | 0.64803139 |
| multipass | TreeFam-A | TPR | 0.66235705 | PPV | 0.69674131 |
| multipass | FAS | NR_ORTHOLOGS | 32012674.0 | FAS | 0.7737068425899689 |
| multipass_refined | GO | NR_ORTHOLOGS | 144241 | avg Schlicker | 0.47253671 |
| multipass_refined | EC | NR_ORTHOLOGS | 188162 | avg Schlicker | 0.93117452 |
| multipass_refined | VGNC | TPR | 0.83274839140971 | PPV | 0.5579631029366479 |
| multipass_refined | SwissTrees | TPR | 0.7104991 | PPV | 0.64755256 |
| multipass_refined | TreeFam-A | TPR | 0.44791715 | PPV | 0.80963845 |
| multipass_refined | FAS | NR_ORTHOLOGS | 8710340.0 | FAS | 0.7752526000932909 |
| strict_profiles | GO | NR_ORTHOLOGS | 351531 | avg Schlicker | 0.44400634 |
| strict_profiles | EC | NR_ORTHOLOGS | 313912 | avg Schlicker | 0.88441954 |
| strict_profiles | VGNC | TPR | 0.9687473886521267 | PPV | 0.16028370755447408 |
| strict_profiles | SwissTrees | TPR | 0.71709513 | PPV | 0.64508979 |
| strict_profiles | TreeFam-A | TPR | 0.6627704 | PPV | 0.69583142 |
| strict_profiles | FAS | NR_ORTHOLOGS | 31682321.0 | FAS | 0.7715748893898831 |
| strict_profiles_refined | GO | NR_ORTHOLOGS | 144347 | avg Schlicker | 0.47235028 |
| strict_profiles_refined | EC | NR_ORTHOLOGS | 187772 | avg Schlicker | 0.93134728 |
| strict_profiles_refined | VGNC | TPR | 0.8321216679201137 | PPV | 0.5575275740440065 |
| strict_profiles_refined | SwissTrees | TPR | 0.70586947 | PPV | 0.64461096 |
| strict_profiles_refined | TreeFam-A | TPR | 0.44803197 | PPV | 0.80770027 |
| strict_profiles_refined | FAS | NR_ORTHOLOGS | 8710722.0 | FAS | 0.7746333592919895 |

## Limitations

- Development-exposed cluster-derived predictions, not native reconciled ortholog pairs.
- Profile contrasts include downstream graph and singleton-assignment responses.
- Refined means sequence-based post-clustering refinement, not phylogeny.
- Six-metric mean is a project-defined secondary summary, not official F1.
- Native standard errors do not supply paired difference confidence intervals.
- No significance, independent generalization, superiority or controlled timing claim.
