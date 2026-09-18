# SwissTrees Relations Involving Missing Input Identities

The frozen-input alias audit identified14reference numeric protein identities
with no accession among the78input FASTAs. This follow-up counts their
incident relations in all eight retained native SwissTrees scoring outputs.
It checks exact common reference identities/truth, family membership, saved
confusion counts and raw file hashes. No score is recomputed or replaced.

```sh
python benchmark_tools/audit_swiss_missing_input_relations.py \
  --counts benchmark_tools/results/qfo_swiss_comparator_counts_20260917.json \
  --aliases benchmark_tools/results/swiss_sequence_alias_audit_20260917.json \
  --output benchmark_tools/results/swiss_missing_input_relations_20260917.json
```

Report SHA-256:
`a8c1749ed6549e71a9aa94af5bb2fa04be5a36b11b3539219250a531be65e4b3`.

## Observed Counts

All eight methods have exactly the same affected labels:181FN,396TN,
zeroTP and zeroFP. The577affected relations are among10,765total reference
relations. Fifteen have both endpoints missing and are counted once.

| Family | Affected FN | Affected TN | Affected Relations |
| --- | ---: | ---: | ---: |
| APP | 22 | 29 | 51 |
| ASTER | 31 | 0 | 31 |
| HOX | 24 | 241 | 265 |
| NOX | 48 | 96 | 144 |
| POP | 10 | 23 | 33 |
| VATB | 46 | 7 | 53 |
| Total | 181 | 396 | 577 |

The other12families have no affected relations. The complete18-family,
eight-method count matrix and original raw-file identities are in the JSON.
Counts are native raw one-direction relations, before the scorer's half-count
and plus-one prior. They are not independent bootstrap observations.

## Interpretation

This measures a common false-negative component in the retained comparison
outputs. It does not explain method-specific predictions on these relations:
every method has the same labels. Removing a shared subset could nevertheless
change aggregate differences because the benchmark averages family precision
and recall and then takes their harmonic mean. Ranking invariance is not
established, and no counterfactual score or improved-input performance is claimed.

Input absence is directly verified for the frozen factorial FASTAs, not
independently reconstructed historical execution inputs. The common labels
are consistent with shared missing inputs but do not prove historical input
parity or the cause of the resource mismatch. Complete those provenance
checks before attributing every shared FN to an identical preprocessing issue.

Retain official full-reference endpoints and their denominators. A separately
prespecified available-input sensitivity analysis, if performed, must remain
secondary and must not silently replace primary evidence. This audit changes
neither inputs nor scoring. Five focused tests and nine existing raw-reader
tests pass, covering each label, double-missing endpoints, empty incidence,
truth/identity changes and duplicate pairs.
