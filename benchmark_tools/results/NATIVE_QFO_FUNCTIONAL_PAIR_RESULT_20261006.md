# Native QfO Functional-Pair Composition

The new native-only diagnostic completes after preparation commit `fe6fe1af`,
without new inference, annotation scoring, admission or bootstrap. The
[machine-readable result](native_qfo_functional_pair_composition_20261006_v1.json)
is 12,376 bytes, SHA256
`1807b3269518590eb339637cf618f07c99abdd7e547e268fc0fea932b4ddf853`.
It exactly replays the admitted snapshot in original Python 3.10 and binds
six raw tables through original admissions and assessment output records.

| Endpoint | R-Off Scored Pairs | R-On Scored Pairs | Shared Scored Pairs | R-Off Native Mean | R-On Native Mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| GO | 145,142 | 78,607 | 78,607 | 0.47211905 | 0.49025964 |
| EC | 186,098 | 116,929 | 116,929 | 0.93211406 | 0.96770222 |
| FAS realized sample | 38,205 | 252,451 | 1,007 | 0.7749447305326528 | 0.7850068177983054 |

Scored GO/EC counts are challenge-eligible relations, not all predictions.
FAS rows are realized scored samples, not the eligible-population denominator
(9,009,082 R-off and 5,113,820 R-on). Neither denominator is interchangeable
with native inference input coverage or an orthology F1 statistic.

## GO And EC

Every R-on scored pair is present in the R-off scored set, with identical
six-decimal scores on all shared pairs. R-off has 66,535 GO and 69,169 EC
scored pairs absent from R-on; their respective rounded-score means are
0.4506870677 and 0.8719528955, lower than the corresponding shared means.
Serialized endpoint differences are +0.0181405879 GO and +0.0355881593 EC.
At retained precision, the increases are entirely a scored-pair membership
and denominator difference, not changes in similarity values on common pairs.
This does not prove identical unavailable full-precision values, annotation
completeness, correct ortholog selection or a biological causal mechanism.

The original-denominator decomposition is preserved in the JSON. Removing
low-scoring relations can raise a functional mean without establishing
improved orthology accuracy. Do not substitute intersection-only means,
which condition away the selection difference and change the endpoint.

## FAS

Only 1,007 realized sampled pairs overlap: 0.3989% of R-on and 2.6358% of
R-off. Their serialized values are identical, with conditional shared-pair
mean difference zero. The original sample mean difference remains
0.010062087265652608; it is not replaced by the shared-only statistic.
The unseeded, method-specific precomputed/missing-score sample mixture and
omitted missing results prevent this small shared subset from representing
either full eligible population. Sample overlap alone cannot determine bias,
recover inclusion probabilities or establish paired confidence intervals.

| Raw Table | Distinct Proteins | Proteins In Multiple Pairs | Maximum Protein Degree |
| --- | ---: | ---: | ---: |
| R-off GO | 47,278 | 39,408 | 44 |
| R-on GO | 43,771 | 32,725 | 38 |
| R-off EC | 26,569 | 24,938 | 63 |
| R-on EC | 25,601 | 22,927 | 54 |
| R-off FAS sample | 66,678 | 8,506 | 6 |
| R-on FAS sample | 257,077 | 126,521 | 17 |

These are dependency diagnostics, not validated independent resampling units.
Native GO/EC Student-t half-widths and FAS pair-IID SEM do not establish
method-difference intervals. Other-endpoint paired uncertainty remains open.

## Independent Readback And Scope

An [independent SQLite join](native_qfo_functional_pair_sql_readback_20261006.json)
parses the six raw tables without the primary parsers/comparison helper.
All 817,432 rows, exact integer GO/EC totals, shared counts and score equality,
FAS overlap and floating-point means match. The receipt is 4,126 bytes, SHA256
`db4f3b6771a4627f5ae039c489a63f9fe6e35c2ad47495d971fdfdc81ee63e24`.
126 joined tests pass in 2.02 seconds, including 13 new independent-join
cases and the 21 primary native-composition cases. Fixtures do not replace
actual provenance replay; the subsequent real invocation checks retained data.

The [GNU-time receipt](native_qfo_functional_pair_diagnostic_resources_20261006.txt)
records only this read-only diagnostic: 3.41 seconds elapsed and maximum RSS
181,976 KiB (about 177.7 MiB), exit zero. It is not native inference timing,
a scaling repeat or a tool-speed comparison. This is a shared-Threadripper
observation under competing work, with unknown tool-dependent effects.
Original failed 22437 timing remains null/ineligible; no timing repair follows.

Reproduce the primary diagnostic using the
[prospective protocol](NATIVE_QFO_FUNCTIONAL_PAIR_PROTOCOL_20261006.md), then
read it back to a fresh destination:

```bash
benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.readback_native_qfo_functional_pairs \
  --composition benchmark_tools/results/native_qfo_functional_pair_composition_20261006_v1.json \
  --composition-sha256 1807b3269518590eb339637cf618f07c99abdd7e547e268fc0fea932b4ddf853 \
  --output /tmp/native_qfo_functional_pair_new_readback.json
```

This is a P0/C0 native ablation, not the selected-default all-tool comparison;
initial HMM search remains on. Historical all-tool diagnostics, endpoint
scores, prior archive/PDF payloads and frozen scientific helpers are unchanged.
The full publication goal remains active and readiness unproven.

After manuscript/checklist/guide integration, 237 joined tests pass in 9.49
seconds, including actual artifact/source/receipt and documentation contracts.
This does not repeat the scientific admission, historical all-tool overlap
panel, large resource census, native inference or full publication audit.
