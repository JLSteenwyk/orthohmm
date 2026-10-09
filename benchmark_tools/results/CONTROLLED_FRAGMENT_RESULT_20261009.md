# Controlled Fragment Observation Results

The single sequential job24087 completed with exit0:0 in7min29s on the shared
Threadripper. All30 new native inference processes and10 diagnostic checkpoints
passed native admission. Retained controls were not rerun. No automatic retry,
failed native identity, missing method outcome or parameter change occurred.

The prespecified condition retained the center60% of floor20% of extant genes
selected solely by seeded ID hashes:1636 observations among8200 genes over
ten independent simulation master seeds. All IDs, species owners and exact
evolutionary-truth bytes were preserved. See the frozen protocol/preparation
addendum for selection, exposure and the prospective runtime inventory change.

## Primary Means

These are arithmetic means of ten seed-level cross-species pair micro scores,
not a pooled-pair statistic. Each method has10/10 admitted seeds. Percent units.

| Method | Baseline F1 | Fragment F1 | Fragment Precision | Fragment Recall |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 97.077 | 96.977 | 94.904 | 99.186 |
| OrthoHMM satellite_v2 phylogeny | 99.458 | 99.336 | 99.894 | 98.806 |
| OrthoFinder3.1.5 full | 99.942 | 99.908 | 99.926 | 99.891 |
| OrthoFinder MCL checkpoint, diagnostic | 97.794 | 97.426 | 95.088 | 99.939 |

The MCL checkpoint is taken from the full run, not a separately executed
sequence-only pipeline or independently timed comparator.

## Paired Uncertainty

All15 planned F1/precision/recall comparisons used all ten paired seeds and
20,000 whole-seed PCG64 draws, seed20261011. F1 differences below are percentage
points. Fixed15-endpoint Bonferroni percentile bounds are conditional
approximations, not exact simultaneous guarantees or independent biological
generalization. All precision/recall intervals are retained in the full table.

| F1 Comparison | Difference | Adjusted Lower | Adjusted Upper |
| --- | ---: | ---: | ---: |
| Fragment minus baseline, high sensitivity | -0.100 | -0.388 | 0.000 |
| Fragment minus baseline, satellite_v2 | -0.122 | -0.484 | 0.089 |
| Fragment minus baseline, full OrthoFinder | -0.034 | -0.124 | 0.000 |
| Fragment high sensitivity minus full OrthoFinder | -2.931 | -5.682 | -1.086 |
| Fragment satellite_v2 minus full OrthoFinder | -0.572 | -1.981 | -0.009 |

Both OrthoHMM configurations have lower F1 than full OrthoFinder in this fixed
synthetic condition. The adjusted F1 interval is negative for both comparisons.
The small mean fragment-minus-baseline decreases have intervals containing
zero; that does not establish equivalence or absence of fragment sensitivity.
Do not turn this test into a tuning campaign or a general superiority claim.

## Endpoint Strata

The zero/one/two-fragment endpoint bins partition ALL truth and prediction
pairs, including between-origin false positives. The following are summed
count changes across ten seeds, not the primary mean or bootstrap statistic.

| Method | Fragment Endpoints | TP Change | FP Change | FN Change |
| --- | ---: | ---: | ---: | ---: |
| High sensitivity | 0 | -3 | 0 | 3 |
| High sensitivity | 1 | -42 | 0 | 42 |
| High sensitivity | 2 | -12 | 0 | 12 |
| Satellite_v2 | 0 | -3 | 5 | 3 |
| Satellite_v2 | 1 | -47 | 1 | 47 |
| Satellite_v2 | 2 | -13 | 0 | 13 |
| Full OrthoFinder | 0 | 2 | 4 | -2 |
| Full OrthoFinder | 1 | -18 | 0 | 18 |
| Full OrthoFinder | 2 | 0 | 0 | 0 |
| MCL checkpoint | 0 | 0 | 153 | 0 |
| MCL checkpoint | 1 | -18 | 63 | 18 |
| MCL checkpoint | 2 | 0 | 6 | 0 |

Untruncated-endpoint pairs can also change through grouping/reconciliation;
the observed differences are not restricted mechanically to incident pairs.
These count patterns alone do not identify the causal pipeline stage.

All20 fragment OrthoHMM native metrics record98-100 constructed expansion
profiles and zero added profile-expansion edges. This confirms execution but
does not demonstrate an accuracy contribution from that expansion stage.
Initial search HMMs remain on; no HMM-off causal comparison is introduced here.
Actual resolved search parallelism is one worker with four threads and a
CPU budget of4, despite the scheduler's16-CPU cap.

## Validation And Scope

Tested source2125e688 was pushed before the one terminal assembly. Results
contain80 baseline/fragment score records,240 stratum records and15 comparisons.
Independent readback used direct native pair/group expansion, Fraction-based
counts/ratios, and separately sorted linear bootstrap interpolation. Only the
existing audited MCL syntax parser was shared. It checked all records and TSV
values without importing the new producer, transform or scoring modules.

Machine-readable report:
`controlled_fragment_results_20261009_v1/report.json`, SHA256
`b1480e8eb2ace2aed743f5fc4e515e27ac4ba68b48e644b33989eda3cf8a7ab7`.
The adjacent scores.tsv/strata.tsv/comparisons.tsv preserve full raw values.
Actual checks, failed historical-inventory preflight, prospective amendment,
job accounting and independent reader command are recorded in
`controlled_fragment_execution_20261009_v1.json`.

This is a development-exposed synthetic observation diagnostic, not natural
fragment truth, domain/indel evolution or independent biological confirmation.
The current OrthoHMM inventory explicitly differs from historical controls;
required scientific versions, sources and native binaries are unchanged,
but historical output equivalence has not been proven. Shared-host resource
observations are potentially tool-dependently confounded; no isolated ranking
or paired comparison to cached baseline timings is justified.

The full publication goal remains incomplete. Next integrate this verified
result, uncertainty, all endpoint strata and its limits into the manuscript
and claim register; preserve the prior frozen manuscript and evidence.
