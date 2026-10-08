# Composed-Native Manuscript Integration

## Prospective Scope

`prepare_composed_native_main_text.py` prepares a separate v5 source and a
native claim-to-evidence table from actual admitted evidence. No production
v5 is generated while native12inference/scoring/admission is unfinished.
Synthetic test documents are not publication results. Preserve every frozen
v4 byte, source, figure, review and archive; no old receipt establishes review
of a new source.

The frozen parent digest is
`0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a`.
The actual retained native11failure-addendum digest is
`4f1a5c536af8ce992d46866374b3ab44ee3158077316d69ba2e9f5fb34cee203`.
Both are mandatory. Other inputs are explicit actual path/digest pairs for
the composed five-cell snapshot, newly computed native-count bootstrap,
production figure manifest and independent figure readback. Source identities
and cross-artifact snapshot/interval/manifest/failure bindings must match.

## Changes And Invariants

- Update only the manuscript header, the affected Abstract sentences, the
  native-QfO section, one Methods paragraph and affected availability links.
- Preserve unrelated all-tool, OrthoBench, simulation, generalization,
  biological examples, error strata, Discussion and References content.
- Retain30admitted endpoints, two missing score cells and all seven coverage/
  status outcomes. Coverage from failed native11scoring is not accuracy.
- Retain FAS137/OUT_OF_MEMORY,8CPU/32GiBfailed-scoring allocation, zero admitted
  endpoint scores and unknown scoring peak memory. Successful native12scoring
  uses its separate8CPU/128GiBallocation, not the matched inference panel.
- Report the six-metric mean only as a project-defined secondary summary;
  validate its arithmetic. GO/EC/FAS are not F1.
- Report actual native-count intervals and100000new shared-family draws with
  seed20260922and all42adjusted endpoints. The manuscript generator itself
  makes zero draws. Missing native cells do not supply cached predictions.
- Derive the numbers of strictly positive, strictly negative and zero-including
  adjusted F1 intervals from the supplied data, including exact-zero bounds.
  Do not copy historical zero-crossing or superiority assumptions.
- Preserve the old complete profile-pair localization verbatim, explicitly
  scoped to P0/C0/R1versusP1/C0/R1. It does not diagnose native12by assertion.
- Generate five native claim-to-evidence rows with explicit excluded claims.
  This is an addition, not a replacement for the full-study claim checklist.
- Retain development exposure, family exchangeability, approximate percentile
  coverage, missing other-endpoint uncertainty, family-overlap and shared-host
  timing limitations. No independent confirmation, new default, isolated
  efficiency, completed full factorial or general superiority is claimed.

All output paths must be distinct, fresh and beside the evidence directory so
relative links remain meaningful. Input/source/evidence hashes are checked
before and after generation. Existing outputs are refused before input reads.

## Commands After Actual Production Figure Readback

Use the retained figure Python environment with Matplotlib/NumPy, disable
Python/user-site/dynamic-loader injection, and constrain BLAS/OpenMP to one
thread. Set all paths/digests from the actual produced artifacts, not test data.

```bash
python -B -m benchmark_tools.prepare_composed_native_main_text \
  --evidence-directory benchmark_tools/results \
  --parent benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md \
    0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a \
  --snapshot "$SNAPSHOT" "$SNAPSHOT_SHA" \
  --intervals "$NATIVE_INTERVALS" "$NATIVE_INTERVALS_SHA" \
  --figure "$FIGURE_MANIFEST" "$FIGURE_MANIFEST_SHA" \
  --reader "$FIGURE_READBACK" "$FIGURE_READBACK_SHA" \
  --failure benchmark_tools/results/native11_qfo_scoring_failure_addendum_20261008_v1/report.json \
    4f1a5c536af8ce992d46866374b3ab44ee3158077316d69ba2e9f5fb34cee203 \
  --output "$FRESH_V5_MARKDOWN" --receipt "$FRESH_GENERATION_RECEIPT" \
  --claims-output "$FRESH_NATIVE_CLAIMS_TSV"
```

If final12cannot be admitted, do not generate this five-cell revision. Preserve
the existing four-cell figure and report the actual failure/missing outcome
in a separate honest addendum. No inference or scoring retry is authorized.

## Rendering And Delivery

Inspection found that the existing `render_manuscript_review.py` and
`print_manuscript_review.py` already accept explicit manuscript/bibliography/
fresh output paths. They do not need a compatibility rewrite for v5.
Reuse these unchanged helpers after actual generation, with the retained
IQ-TREE3-updated October bibliography, not an older September selection.
Validate citations and links, print into a fresh directory, and independently
inspect the actual new PDF/HTML and every figure. Do not assume the prior
page count, prior table completeness or prior visual review transfers.
Generate a new direct-review selection/package from the actual new receipts,
without rebuilding historical archives or claiming transitive study closure,
redistribution clearance or an archival DOI.

## Test Meaning

Fixtures explicitly invent the final admission and use real count/bootstrap,
render, asset-readback, generation and CSV kernels. Tests demonstrate targeted
manuscript replacement and byte-equivalent unrelated text sections; positive,
negative and exact-zero interval interpretation; refusal of mixed sources,
cross-artifact bindings, partial readbacks, missing/boolean/nonfinite means,
invented failure peak memory, altered failure allocation, bad anchors and
existing/aliased/outside outputs. No fixture is a production native12result.
