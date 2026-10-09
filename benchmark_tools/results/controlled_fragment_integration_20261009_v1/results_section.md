### Controlled Fragment Accuracy

All 40 fragment method/checkpoint outcomes were admitted across ten seeds.
HS denotes high sensitivity, PHY denotes satellite_v2 with phylogeny, and
OF denotes full OrthoFinder. Values below are percent; coverage is defined
in the controlled observation protocol above.

| Method | Baseline F1 | Fragment F1 | Fragment Precision | Fragment Recall | Coverage |
| --- | ---: | ---: | ---: | ---: | ---: |
| HS | 97.0772 | 96.9773 | 94.9040 | 99.1855 | 99.9635 |
| PHY | 99.4584 | 99.3361 | 99.8940 | 98.8060 | 99.9147 |
| OF | 99.9424 | 99.9080 | 99.9257 | 99.8906 | 99.9633 |
| OF checkpoint | 97.7938 | 97.4260 | 95.0881 | 99.9386 | 99.9878 |

### Controlled Fragment Paired Differences

All 15 planned differences and interval endpoints are shown in percentage
points. F-B denotes fragment minus baseline; F-F denotes the fragment
method minus fragment OF. Each row has ten eligible paired seeds.

| Comparison | Metric | Difference | Nominal 95% Interval | Fixed-15 Adjusted Interval | Seeds |
| --- | --- | ---: | --- | --- | ---: |
| HS F-B | f1 | -0.0999 | [-0.2759, +0.0000] | [-0.3878, +0.0000] | 10 |
| HS F-B | precision | -0.0188 | [-0.0539, +0.0000] | [-0.0740, +0.0000] | 10 |
| HS F-B | recall | -0.1891 | [-0.5196, +0.0000] | [-0.7326, +0.0000] | 10 |
| PHY F-B | f1 | -0.1223 | [-0.3502, +0.0377] | [-0.4835, +0.0893] | 10 |
| PHY F-B | precision | -0.0222 | [-0.0836, +0.0224] | [-0.1126, +0.0376] | 10 |
| PHY F-B | recall | -0.2116 | [-0.6017, +0.0815] | [-0.8371, +0.1805] | 10 |
| OF F-B | f1 | -0.0344 | [-0.0893, +0.0000] | [-0.1237, +0.0000] | 10 |
| OF F-B | precision | -0.0139 | [-0.0347, +0.0000] | [-0.0417, +0.0000] | 10 |
| OF F-B | recall | -0.0547 | [-0.1503, +0.0000] | [-0.2050, +0.0000] | 10 |
| HS-OF F-F | f1 | -2.9307 | [-4.6881, -1.5314] | [-5.6825, -1.0863] | 10 |
| HS-OF F-F | precision | -5.0217 | [-7.2106, -3.0158] | [-8.2208, -2.1726] | 10 |
| HS-OF F-F | recall | -0.7051 | [-2.2931, +0.1571] | [-3.0668, +0.2186] | 10 |
| PHY-OF F-F | f1 | -0.5719 | [-1.3924, -0.0404] | [-1.9815, -0.0087] | 10 |
| PHY-OF F-F | precision | -0.0316 | [-0.1206, +0.0188] | [-0.1655, +0.0236] | 10 |
| PHY-OF F-F | recall | -1.0846 | [-2.5981, -0.0906] | [-3.6622, -0.0244] | 10 |

Both fragment HMM configurations have lower mean F1 than full OF in this
condition; both adjusted F1 difference intervals exclude zero on the
negative side. The three native fragment-minus-baseline F1 intervals
include zero, which establishes neither equivalence nor absence of
fragment sensitivity. No default was changed based on these outcomes.

### Controlled Fragment Endpoint Counts

Counts below sum the ten seeds; they are not the mean-score or interval
statistics. Endpoint count is the number of prospectively flagged genes
in a pair, including the reused untruncated baseline arm.

| Method | Truncated Endpoints | Baseline TP / FP / FN | Fragment TP / FP / FN |
| --- | ---: | --- | --- |
| HS | 0 | 18538 / 1017 / 116 | 18535 / 1017 / 119 |
| HS | 1 | 9169 / 521 / 66 | 9127 / 521 / 108 |
| HS | 2 | 1155 / 54 / 7 | 1143 / 54 / 19 |
| PHY | 0 | 18487 / 14 / 167 | 18484 / 19 / 170 |
| PHY | 1 | 9122 / 10 / 113 | 9075 / 11 / 160 |
| PHY | 2 | 1148 / 0 / 14 | 1135 / 0 / 27 |
| OF | 0 | 18639 / 13 / 15 | 18641 / 17 / 13 |
| OF | 1 | 9234 / 5 / 1 | 9216 / 5 / 19 |
| OF | 2 | 1162 / 0 / 0 | 1162 / 0 / 0 |
| OF checkpoint | 0 | 18654 / 853 / 0 | 18654 / 1006 / 0 |
| OF checkpoint | 1 | 9235 / 430 / 0 | 9217 / 493 / 18 |
| OF checkpoint | 2 | 1162 / 48 / 0 | 1162 / 54 / 0 |

Unflagged-pair counts can change through grouping and reconciliation even
though neither endpoint was truncated. These counts alone do not localize
a causal stage. Representative search/candidate/group/reconciliation tracing
remains separate work. The initial HMM search is on in both HMM methods;
this observation test is not an initial-HMM-off causal comparison.
[All 80 scores, 240 strata and 15 comparisons](controlled_fragment_results_20261009_v1/report.json),
[full-precision stratum table](controlled_fragment_results_20261009_v1/strata.tsv),
[result interpretation](CONTROLLED_FRAGMENT_RESULT_20261009.md).

