### Shared-Host Resource Panel

The snapshot retains 27 of 27 reviewed attempts, 25 with measured native resources and 24 eligible observations. 6 method/size cells have three eligible repeats.

These are shared-host matched-resource observations with 32 native CPUs and a 128-GiB limit per run, not isolated comparative timing. Contention distortion is unknown and potentially method dependent. No background overhead is subtracted.

Tables report median [minimum, maximum] only for cells with three eligible repeats. An unavailable summary is not zero or a failed native inference. Ranges are not confidence intervals. Eligibility counts retain measured exclusions and pre-native aborts rather than selecting the fastest attempts.

#### Native command wall (seconds)

| Method | Proteomes | Eligible/planned | Median [minimum, maximum] |
| --- | --- | --- | --- |
| OrthoHMM high sensitivity | 4 | 2/3 | Unavailable |
| OrthoHMM high sensitivity | 8 | 3/3 | 1238.9352 [1221.4444, 1303.6330] |
| OrthoHMM high sensitivity | 12 | 3/3 | 2720.0774 [2621.4976, 3027.0605] |
| OrthoHMM inferred phylogeny | 4 | 2/3 | Unavailable |
| OrthoHMM inferred phylogeny | 8 | 3/3 | 1999.5508 [1888.8847, 2086.9473] |
| OrthoHMM inferred phylogeny | 12 | 3/3 | 4392.8716 [4042.4262, 4809.2924] |
| OrthoFinder 3.1.5 full | 4 | 3/3 | 441.2890 [420.7400, 474.2999] |
| OrthoFinder 3.1.5 full | 8 | 3/3 | 1139.9398 [1120.2548, 1179.9703] |
| OrthoFinder 3.1.5 full | 12 | 2/3 | Unavailable |

#### Task-subtree CPU bracket (CPU-seconds)

| Method | Proteomes | Eligible/planned | Median [minimum, maximum] |
| --- | --- | --- | --- |
| OrthoHMM high sensitivity | 4 | 2/3 | Unavailable |
| OrthoHMM high sensitivity | 8 | 3/3 | 34395.7481 [33878.9946, 36610.5696] |
| OrthoHMM high sensitivity | 12 | 3/3 | 78335.8306 [75126.8588, 86415.3192] |
| OrthoHMM inferred phylogeny | 4 | 2/3 | Unavailable |
| OrthoHMM inferred phylogeny | 8 | 3/3 | 56008.3471 [52911.7819, 58240.7234] |
| OrthoHMM inferred phylogeny | 12 | 3/3 | 125264.7375 [116000.9964, 136293.3322] |
| OrthoFinder 3.1.5 full | 4 | 3/3 | 4602.7585 [4455.3737, 5008.8970] |
| OrthoFinder 3.1.5 full | 8 | 3/3 | 17734.4088 [17232.2843, 18787.3118] |
| OrthoFinder 3.1.5 full | 12 | 2/3 | Unavailable |

#### Native-step lifetime peak (GiB)

| Method | Proteomes | Eligible/planned | Median [minimum, maximum] |
| --- | --- | --- | --- |
| OrthoHMM high sensitivity | 4 | 2/3 | Unavailable |
| OrthoHMM high sensitivity | 8 | 3/3 | 7.724 [7.715, 7.932] |
| OrthoHMM high sensitivity | 12 | 3/3 | 10.555 [10.533, 10.730] |
| OrthoHMM inferred phylogeny | 4 | 2/3 | Unavailable |
| OrthoHMM inferred phylogeny | 8 | 3/3 | 7.681 [7.677, 7.702] |
| OrthoHMM inferred phylogeny | 12 | 3/3 | 10.511 [10.503, 10.603] |
| OrthoFinder 3.1.5 full | 4 | 3/3 | 6.775 [6.769, 6.928] |
| OrthoFinder 3.1.5 full | 8 | 3/3 | 10.226 [10.215, 10.233] |
| OrthoFinder 3.1.5 full | 12 | 2/3 | Unavailable |

Excluded attempt indices: 0, 17, 20. Pre-native abort indices: 17, 20; these have no native resource endpoints, not zero measurements. Unreviewed indices: none.

Across measured attempts, maximum observed foreign CPU demand ranges from 42.1424 to 72.8708 core equivalents. These are process-interval observations, not reservations or estimates of causal slowdown.

CPU includes the native task-subtree wrapper bracket. Peak memory includes the native-step launcher and is not pure algorithm RSS. Preparation, conversion and scoring are outside the native command timer. This single nested taxon series co-varies proteome count and taxon composition, not taxon-invariant scaling. Native output validation does not provide new prediction-accuracy evidence. Historical DGX and earlier shared-host timings are not pooled.
