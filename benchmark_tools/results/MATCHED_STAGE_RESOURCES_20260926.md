# Descriptive Matched-Control Stage Resources

[Machine-readable records](matched_stage_resources_20260926.json) retain
each stage's checksummed GNU-time log, elapsed/user/system time, reported
maximum process RSS and exit code. All 35 reporting datasets have both search
and graph observations. DIAMOND database preparation is reported separately.
The strict verbose-log adapter reuses the existing GNU-time numeric/semantics
validator; all 19 focused tests pass.

Median observations across reporting datasets:

| Arm / phase | Elapsed seconds | CPU seconds | Maximum process RSS (KiB) |
|---|---:|---:|---:|
| HMM search | 5.18 | 19.43 | 119944 |
| HMM graph | 0.24 | 0.23 | 41472 |
| DIAMOND database preparation | 0.00 | 0.01 | 13824 |
| DIAMOND broad search | 1.88 | 5.68 | 21504 |
| DIAMOND graph | 0.23 | 0.22 | 39936 |

The JSON includes minima/maxima and individual stage measurements. Phase
elapsed/CPU values sum sequential native commands within each dataset;
RSS takes the maximum reported process value across those commands, not
their sum. GNU time rounds these small measurements: 0.00 does not mean
database preparation is free. RSS is not concurrent process-tree or cgroup
memory. Wrapper, hashing and copying overhead are excluded.

These are shared-host, incremental observations, **not controlled timing**.
DIAMOND searched at E=1 and was subsequently filtered to E<=1e-40; search
effort was not matched. Neither a comparative speedup nor memory-efficiency
claim is admitted. Offline numeric-export preparation and scoring were not
separately timed, so no complete end-to-end resource total is supplied.
The requirement for dedicated-machine, matched-resource scaling remains open.

These observations complement the accuracy/coverage results without changing
their conclusions or scientific defaults. All raw inference and scoring
artifacts are retained; no native job was repeated for this summary.
