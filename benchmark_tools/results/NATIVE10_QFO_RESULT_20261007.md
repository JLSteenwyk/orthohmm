# Native10 QfO Result

This chronological result supplements the frozen three-cell manuscript and
figures; it does not rewrite their historical counts. The
[four-cell snapshot](native_qfo_scientific_scores_20261007_v3/scores.md)
now contains four admitted accuracy rows among seven planned QfO identities.
Identity9 remains failed; identities11/12 have no admitted score. Missing
results are not zeros. P1 is downstream profile refinement, C0 disables
candidate expansion, and R1 enables reconciliation; initial HMM search
remains enabled throughout the factorial.

## Independent Admission

Native23902 and distinct full review23973 completed0:0. The original failed
review23910 remains retained with unestablished cause; neither diagnostic
replaces it or changes its bytes. Conversion23977 preserved all5,115,410 native
phylogenetically inferred pairs with zero mapping loss. Assessment23978
completed0:0 in35:13, and admission23984 completed0:0 in3:00, ended15:20:18.
The frozen independent admitter checked all six endpoints, full trace,
output/source bindings and FAS sample membership;113,273 checked records
are retained in the external admission receipt:
`benchmarks/results/allocated_native_qfo_admission_v1/p1_c0_r1/results.json`,
SHA256 `9fa1f24a45e5ea604cc593395acaac9a40ef1189323c6682d4a007a420daffd9`.
Its status is `allocated_native_factorial_qfo_assessment_admitted` and
`accuracy_admitted=true`; publication readiness remains false.

| Endpoint | Score | Precision | Recall | Statistic |
| --- | ---: | ---: | ---: | --- |
| VGNC | 0.896981 | 0.999384 | 0.813612 | F1 |
| SwissTrees | 0.786372 | 0.949051 | 0.671303 | F1 |
| TreeFam-A | 0.601468 | 0.957097 | 0.438525 | F1 |
| GO | 0.490409 | Not applicable | Not applicable | Average Schlicker similarity |
| EC | 0.967582 | Not applicable | Not applicable | Average Schlicker similarity |
| FAS | 0.785014 | Not applicable | Not applicable | Feature architecture similarity |

The six-metric mean is0.7546377200341031, a project-defined secondary
summary, not official QfO F1. GO/EC assessed117,068/78,843 relations in EC/GO
order; those endpoint-specific populations are not submitted-pair coverage.
542,336 of984,137 full input accessions occur in at least one predicted
relation (55.1078%); relation coverage is not reference recall or accuracy.

FAS retains257,182 sample pairs/259,532 proteins, sample mean
0.7850143560591332 and native pair-IID SEM0.0003307313856106987. The native
reported eligible-pair count is5,115,410 and sample fraction0.050275930961545603.
128,887 proteins appear in multiple sampled pairs; maximum sample degree18.
Membership and arithmetic were independently checked. The frozen protocol
uses unseeded native sampling, up to9,000 newly computed missing-pair scores,
and omits failed/uncomputed missing scores without zero imputation. Its
eligible population is not guaranteed to be every prediction. The SEM is
not an independent-family or paired-method confidence interval, and verified
sample arithmetic does not prove representativeness.

## Resources And Limits

The retained terminal review supplies53,995.367466727 wall-seconds,
1,585,693.069194 CPU-seconds and19,355,951,104 peak-memory bytes. Wall scope
is the native command monotonic interval; CPU is the native-task subtree
bracket including its wrapper; peak is native-step lifetime including launcher.
Maximum retained foreign average CPU demand is110.64280532815579 cores.
Conversion, scoring, admission and reporting costs are separate. Assessment
elapsed35:13, FAS task30m47s, admission3:00; the direct-metadata exporter
exited0 in2.49s with517,452KiB maximum RSS and zero swaps. These are not new
inference repeats or isolated efficiency measurements. Scientific timing
admission remains false, as specified by the retained contract.

Timing measurements were collected on a shared Threadripper while other
analyses were running. Competition for CPU, memory bandwidth and I/O may
have affected elapsed times, with an unknown and potentially tool-dependent
impact. These are observed shared-host timings, not estimates of isolated
performance.

QfO remains development-exposed. A fourth cell does not complete the
interaction matrix, establish independent confirmation or general superiority
over full OrthoFinder, justify default tuning, or supply unresolved paired
uncertainty for other endpoints/the secondary mean. Native resolved pairs
and R-off group-derived clique pairs have different prediction semantics.
The next genuinely unrun identities remain11(P1C1R0), then12(P1C1R1),
sequentially under the frozen plan and reviewed history, with failure9 retained.

## Reporting Provenance

The existing committed exporter was executed once in scientific Python3.10
without source changes or raw admission replay. The new
[machine-readable snapshot](native_qfo_scientific_scores_20261007_v3/report.json)
has SHA256 `7aab1cbd31fb650df42e6e80a14ee0167b35a810f08202768c89c514755211a2`.
[TSV](native_qfo_scientific_scores_20261007_v3/scores.tsv) SHA256
`aeff48737c87a3e61ef006a17e6de76fffb8400bf84c2d572f5dcb1b871d81da`;
Markdown SHA256 `fccbdbe7f20ff0da70d75a066668a010c008d7966939a1df870d857a01d1c6ed`.
A separate stdlib JSON/CSV readback confirmed all six new score formulae,
mean, denominator, seven-row inventory and missing cells; every prior row
is exactly unchanged from the preceding snapshot. This is direct reporting
validation, not another scientific admission or full-study runtime closure.

Reproduction uses fresh output only; retained outputs must not be overwritten:

```bash
env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.export_allocated_native_qfo_scientific_scores \
  --plan benchmark_tools/results/native_factorial_receipt_amendment_20261004/plan.json \
  --plan-sha256 6c87babcbb5581830e0b9e7b9bf9aaba30a85bde1c4ab465e561017e67e9c89b \
  --admission benchmarks/results/full_native_qfo_admission_v1/p0_c0_r0/results.json \
    1e5b7ba7824237080286a491ef2bf7c2566627bf614c7e667551fa6b566f6411 \
  --admission benchmarks/results/full_native_qfo_admission_v1/p0_c1_r0/results.json \
    7e830778018ae78c32a6c01dee0fa85a413c18320a0d34487244e00971827c60 \
  --recovered-admission benchmarks/results/measurement_failed_native_qfo_admission_v1/p0_c0_r1/results.json \
    4c3a17a76eff5043c8b40d0c9f1e8ead6c7ed32d4e248dbd485988c33537ae3b \
  --allocated-admission benchmarks/results/allocated_native_qfo_admission_v1/p1_c0_r1/results.json \
    9fa1f24a45e5ea604cc593395acaac9a40ef1189323c6682d4a007a420daffd9 \
  --output benchmark_tools/results/native_qfo_scientific_scores_20261007_v3
```
