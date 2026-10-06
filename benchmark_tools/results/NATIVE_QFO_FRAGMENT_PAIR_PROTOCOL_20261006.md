# Native SwissTrees Fragment-Annotation Pair Audit

Specify before calculating fragment-stratified native pair changes. The native
overall scores, complete 334 TP removals/1,689 FP removals and historical
annotation panel have already been inspected. This is retrospective descriptive
error analysis, not a new inferential or independent confirmation endpoint.
Keep the frozen method and all original jobs unchanged.

## Frozen Inputs

- Native transitions `native_qfo_swiss_pair_transitions_20261006_v2.json`, SHA256
  `fe003f2cbc4285ea56cd80a92b71703b5c9e135ad1c8f244504dcfeab0911b46`.
- Historical annotation admission `swiss_historical_fragment_admission_22117.json`,
  SHA256 `a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8`.
- Original membership/features `corrected_swiss_sequence_strata_20260918.json`,
  SHA256 `912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1`.

Use the completed native P0C0 R0/R1 raw SwissTrees files only. Match all 18
families/563 exact accession members to the historical admission and original
recorded FASTA identities. Recheck each of 563 selected `entry.txt` files against
its original admission hash, parse with Bio.SwissProt, and reproduce accession,
taxon, sequence digest/version, entry version and admitted annotation features.
Reuse the prior history-selection admission; do not reacquire histories or
rerun the complete historical annotation admission/whole-proteome extraction.
Report unchecked history/acquisition records as inherited evidence.

## Pair Views And Output

Retain the original sequence-matched historical view and baseline-only
sensitivity view. In baseline-only, the 14 later-sequence-version records count
as missing. A protein is annotation-positive if its admitted `fragment_flag`
or nonempty `incomplete_sequence_features` is positive. An unflagged record is
not proven complete. For each pair use these mutually exclusive bins:

1. `annotation_positive`: at least one positive endpoint.
2. `missing_without_positive`: no positive endpoint and at least one missing.
3. `all_matched_unflagged`: both matched and no positive endpoint.

These are endpoint-level pair bins, not the earlier family-level bins: a
positive family does not make every member or pair annotation-positive.
Positive overrides missing exactly as in the retained historical family rule.

Retain every one of the 10,765 reference relations, unchanged decisions too,
in a sorted ledger with both cell labels and both annotation-view bins. Report
all 16 R0-to-R1 transitions and before/after TP/FP/FN/TN counts per view/bin,
including explicit zero cells/empty bins. Descriptive TP-removal fractions use
R0 TP as denominator, FP-removal fractions R0 FP; zero denominator is null,
not zero. These counts/rates are not native F1, independent-pair uncertainty
or biological validation. Do not manufacture significance or compare rates
as independent binomial trials.

Use a separate reader, importing no primary exporter, to reparse raw relations
and annotation entries, construct an independent SQLite join/grouping, and
verify the full ledger, every transition/marginal and rate. Annotation parses
share Bio.SwissProt: disclose that dependency rather than claim entirely
independent parsers. No new bootstrap draws or transferred historical CIs.

## Interpretation And Delivery

Test the descriptive explanation that true-positive removals principally
involve annotation-positive endpoints by retaining all counts, not selecting
only illustrative cases. No causal attribution to fragments/reconciliation,
experimental completeness claim, tuned exclusion or independent generalization
follows. Original historical records may have annotation errors or temporal
limitations. Failed R1 timing stays ineligible. Both cells keep initial HMM
search on; neither is relabeled as a selected-default configuration.

Commit/push tested sources before actual execution, then generated report,
ledger, independent readback, receipts and bounded manuscript claim links.
Execution costs are shared-host postprocessing only; CPU, memory-bandwidth
and I/O competition have unknown, potentially tool-dependent effects, not
isolated inference-speed evidence. Preserve original protocol, manuscript,
archive and failed/native execution bytes. Full goal remains active.
