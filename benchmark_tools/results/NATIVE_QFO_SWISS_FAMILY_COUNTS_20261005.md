# Native QfO SwissTrees Family Counts Verified

The previous goal turn made concrete progress at pushed `95d8f528`: the
first fresh native QfO scores were admitted and integrated into the extended
manuscript. Re-read the full goal and verify original inference22437 is still
RUNNING, with reviewer22440 dependency-pending. Do not restart either handle
or submit index8 before terminal review and its fresh gates.

## Actual Family Evidence

The [new count audit](native_qfo_swiss_family_counts_20261005_v1.json) reads
the actual inventoried SwissTrees raw file from native22435/conversion22439/
assessment22441/admission22442, not a substituted cached count table. It
replays the supplied-admission score snapshot, checks raw/execution records
were independently inventoried, checks reference labels/members/coverage,
and reconstructs the native family statistics and aggregate. Twenty-two
direct input records are rechecked, separately from the full raw admission.

All 18 native P0C0R0 family records, including raw TP/FP/FN/TN counts,
represented genes and prior-adjusted statistics, exactly match the retained
corrected factorial's corresponding cell. The aggregate also matches:

| Statistic | Full-Precision Count Reconstruction | Native Decimal Endpoint |
| --- | ---: | ---: |
| Precision | 0.6439400942516949 | 0.64394009 |
| Recall | 0.741266050122766 | 0.74126605 |
| F1 of macro precision and recall | 0.689183963152562 | 0.6891839606644314 |

The native scorer uses one-direction raw confusion counts divided by two,
then a prior of one for TP/FP/FN. The statistic is the harmonic mean of mean
family precision and mean family recall, not the mean of family F1. Retain
the native decimal endpoint unchanged; its tiny difference from the
full-precision reconstruction is not a new score or a biological effect.

The [independent readback](native_qfo_swiss_family_readback_20261005_v1.json)
uses a separate CSV/gzip parser and decimal-independent rational arithmetic.
It checks all 10,765 scored reference-pair labels against the retained
corrected cell: identities, TP/FP/FN/TN labels and represented members are
identical. Thus this bounded SwissTrees decision universe matches, not just
the aggregate. This establishes neither whole-partition equivalence nor
equality for any other QfO endpoint or all 9,009,082 submitted relations.

## Reproduction And Tests

Execute the new audit once against the actual committed score snapshot:

```bash
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.audit_native_qfo_swiss_counts \
  --snapshot benchmark_tools/results/native_qfo_scores_20261005_v1/report.json \
  --snapshot-sha256 6bca0b771af5604f76dfe7796fa1590fdfda71e78f7599c28cd6f1555d6a1c13 \
  --retained-counts benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_counts.json \
  --output benchmark_tools/results/native_qfo_swiss_family_counts_20261005_v1.json
```

The original output already exists and the CLI refuses overwrite. For a
genuine reproduction use a fresh destination, not a duplicate native,
converter or assessment job. The command reads existing outputs only.
Actual execution exits0 with one cell and18 family records. A separate
prepared Python3.10 CLI-help probe also passes; this is not fresh-install
or independent biological certification.

All [143 joined tests](native_qfo_swiss_counts_joined_tests_20261005.xml)
pass in2.62s, zero failures/errors/skips, including21 new audit cases. They
exercise actual native reporting bindings, full raw family reconstruction,
changed valid counts retained without substitution, altered truth/member
labels, malformed/duplicated raw rows, missing native metrics, unchecked
execution/raw files, ambiguous raw inventory, snapshot tampering, absent
admissions and no-overwrite behavior. The
[initial receipt](native_qfo_swiss_counts_tests_20261005.xml) retains20 passes
and one failed fixture assertion: the first alternate synthetic cell had
different decisions but identical counts. Select the genuinely count-changing
fixture; no production code or original result is changed to hide failure.

After manuscript/claim/guide integration, the
[final 155-test suite](native_qfo_swiss_integration_tests_20261005.xml) passes
in3.84s, zero failures/errors/skips, adding the12 existing documentation
contracts. Preserve the earlier143-test receipt; neither suite supplies a
biological confirmation or a missing native contrast.

Independent readback verifies the final test counts, source and direct input
pins, all920 unchanged frozen helper files and actual original live PID265470,
creation time1791220584.63/affinity0..31. The new script/test are outside
frozen source pins. No inference/scoring/FAS sampling/resource replay,
bootstrap draw, changed scientific setting or unrelated-work intervention.

## Limits And Next Action

This count audit attaches no historical confidence interval and generates
zero new bootstrap draws. Only one fresh native cell is admitted; six others
remain unavailable, not zero. Full native conditional comparisons require
their own admitted counts. Initial HMM search remains enabled in P-off.
R-on versus R-off changes resolved versus group-clique prediction semantics.
The 18 SwissTrees families remain development-exposed with possible
exchangeability limitations. Other QfO endpoints and the secondary mean
still need their own defensible uncertainty methods. No superiority,
equivalence, independent confirmation or full publication readiness follows.

Observe original22437/22440 next. Preserve the completed first native score,
older cached factorial, failure history and rc4/PDF snapshots. Broader
provenance, uncertainty and final-release requirements remain open. Shared
Threadripper contention has unknown, potentially tool-dependent timing
effects; this count audit is separate from inference timing. No quiet-window
or dedicated-host blocker is introduced. The full goal remains active.
