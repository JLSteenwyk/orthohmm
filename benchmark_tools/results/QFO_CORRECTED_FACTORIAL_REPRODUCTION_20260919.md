# Corrected Factorial Statistical Reproduction

The complete corrected-release SwissTrees factorial statistics reproduced
exactly outside the checkout from committed family counts. Source revision
`0864f9b` exports 15 Git blobs: runner, unchanged shared engine and imports,
counts, admitted result, both frozen protocols, license and analysis
requirements/lock files. Mutable worktree bytes are not scientific inputs.

Execution used `/tmp/orthohmm-corrected-factorial-Mbf9oa/reproduction`, the
existing isolated analysis environment (Python 3.10.13, NumPy 2.2.6), Python
`-I -B`, and one numerical thread. Inherited Python and loader path overrides
were removed. Exact installed distribution versions and exported hashes are
retained in the report. All exported files were rechecked after execution.

The unchanged numerical engine used 100,000 shared family draws, seed
20260922, all 18 families and Bonferroni adjustment over 42 endpoints.
All 11 numerical fields match exactly, including eight cell estimates,
14 contrasts, F1/precision/recall differences, both interval types, family
effects and positive/tie/negative counts. Only the input schema tag is
adapted, as in the admitted corrected wrapper; no counts or arithmetic change.
The output retains the generic engine status and does not repeat source
admission or promote publication readiness.

- [Reproduction report](qfo_corrected_factorial_reproduction_20260919.json),
  SHA-256 `cbdee0c6d499a2a4a6c19edab68ca715d8830b37a718d3ed7484febe02967710`.
- [Admitted scientific results and limitations](QFO_CORRECTED_FACTORIAL_COMPLETE_20260919.md).
- Eighteen focused tests passed, including every numerical field, wrong
  release, admission promotion, changed engine and existing-output rejection.
- The isolated worker exited zero; the report records exact numerical match.

Reproduce with fresh destinations:

```sh
python benchmark_tools/reproduce_corrected_factorial.py \
  --repo . --revision 0864f9b \
  --python benchmarks/work/swiss_analysis_env_20260917/bin/python \
  --output /tmp/corrected-factorial-reproduction \
  --report /tmp/corrected-factorial-reproduction.json
```

The environment must already exist; this run records it but does not reinstall
dependencies or authenticate OS libraries. This is same-host statistical
reproduction, not independent biological confirmation, cross-platform testing,
native inference, pair conversion, official scoring or raw-reference audit.
Historical absolute paths remain provenance and are not accessed as inputs.
Development exposure, limited family exchangeability and all uncertainty
limitations remain unchanged. No new default or superiority claim follows.
