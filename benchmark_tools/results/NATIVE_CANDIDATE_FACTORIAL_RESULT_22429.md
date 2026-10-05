# Native Candidate Order/Score Factorial Result

## Prespecified Execution

Commit6a59e1d0 fixes the [five-arm diagnostic protocol](../NATIVE_CANDIDATE_FACTORIAL_20261004.md).
The first matching-runtime preparation fails before creating any plan or
scientific attempt: an existing data-loader import transitively requires
Numba, absent from the review environment. Commit85c7f41d replaces this
unused search import with checksum-verified, schema-checked local cache
loading. All48 joined tests and the actual Python3.10 CLI import pass before
the real plan is produced. No frozen runtime or scientific source changes.

The actual [plan](native_candidate_factorial_22429/plan.json), SHA
`73edbaee2c1602880ccc0bf5a7ac7d45dcdef007919ee71b819b497324ef3144`,
binds the original seed bytes/order, all five arms/order, satellite_v2
parameters, retained score and both whole partitions, trusted historical
pickle, every native checkpoint file, Python3.10.13/NumPy2.2.6 and source pins.
Available RAM662,731,202,560bytes before launch; only this diagnostic uses
CPU32, leaving native22430's affinity0..31 unchanged. The
[startup record](native_candidate_factorial_22429/started.json) verifies
PID1675877/affinity32 and actual argv/runtime. No Slurm job is launched.

The single attempt exits0 with both complete historical/fresh controls
reproduced. All five arms run without retry, tuning or early termination.
There are18,235,373 unique shared directed nonself hits; the full fresh
control also includes251,135 self hits. All arms cover251,378genes and yield
54,745groups after8,500 merges. The92.67385396000464s duration is descriptive
diagnostic time, not a comparative benchmark measurement.

## Result

The [machine-derived table](native_candidate_factorial_22429/summary.md) and
[exact result](native_candidate_factorial_22429/report.json) show:

| Hit order | Scores | Self hits | Changed genes vs historical | Changed genes vs native |
| --- | --- | --- | ---: | ---: |
| Historical | Historical | Excluded | 0 | 58 |
| Historical | Fresh | Excluded | 0 | 58 |
| Fresh | Historical | Excluded | 58 | 0 |
| Fresh | Fresh | Excluded | 58 | 0 |
| Fresh | Fresh | Included | 58 | 0 |

With seeds, function, parameters and aligned hit values fixed, changing hit
order alone reproduces the entire observed58-gene candidate-partition
difference. Changing score values alone does not alter the historical
partition; changing them under fresh order also does not alter the fresh
partition. Self-hit inclusion has no partition effect in this test.

In the order-only intervention,8,482 semantic accepted merges are common.
For those common merges, hit counts, maximum scores, coverage, species
overlap and source/target numeric IDs are identical. Support changes at
1,725 common merges, by at most1.4210854715202004e-14. Three round0 anchors
select different attachments at cap4; their changed-selection support is
tied or differs by7.105427357601002e-15. The pinned function preserves
within-cluster-pair order with stable sorting before `np.add.reduceat`;
its accumulated support is order-sensitive at this boundary.

## Independent Readback And Limits

The [readback](native_candidate_factorial_readback_22429.json) checks all
five whole accepted-union reconstructions and all ten pair comparisons,
both complete controls and41 direct pins. Its SHA is
`8ad5cd74372c5924c079526ae7301bfce255f981630fca69b2f526626fa434bd`.
The [separate stdlib gene-level union check](native_candidate_factorial_22429/summary.json)
reconstructs all five partitions without importing the replay/partition
helpers, verifies all920 unchanged native helper pins and confirms three
published plan/start/result copies byte-identical to the originals.

The runner pins whole output partitions; merge-trace digests are first
bound by the readback, not falsely attributed to its original per-arm rows.
This is consistency evidence, not independent recalculation of candidate
scores or rejected alternatives. Original raw partitions/traces stay in
the bound local work directory; no large raw files are committed.

All57 joined driver/alignment/replay/readback tests pass, zero failures,
errors or skips,29.21s. The [saved JUnit](native_candidate_factorial_readback_tests_20261004.xml)
includes actual complete-report reproduction, the order-only feature checks,
toy end-to-end execution and negative integrity/failure-handling cases.

The fixed-seed interventions identify an order effect here, not universal
determinism, the original upstream cause of hit order/score bits, exact
fresh seed order, biological accuracy improvement or a future-method
validation. Existing reference accuracy is inherited from the separately
scored unchanged historical/native partitions. Defaults remain frozen;
do not opportunistically round values or choose a matching arm. Shared-host
timing distortion remains unknown and potentially tool-dependent, with no
causal contention claim or timing correction.

Actual22430 remains RUNNING in MAFFT/FastTree inference. Its full terminal
review and separate score still precede the next native identity. Full
publication readiness and remaining scientific requirements are unproven.
