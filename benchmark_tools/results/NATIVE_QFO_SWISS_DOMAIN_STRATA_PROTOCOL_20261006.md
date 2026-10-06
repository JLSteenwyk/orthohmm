# Native SwissTrees Fixed Domain Bins

Specify before computing the native domain-stratified scores. Overall native
scores, family counts and the historical selected-default domain analysis have
already been inspected. This is retrospective development-exposed explanation,
not independent confirmation or a new prespecified inferential test.

## Inputs And Scope

- Native count projection: `native_qfo_swiss_sequence_strata_20261006_v1/report.json`,
  SHA256 `4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3`.
- Independent raw-count reader receipt: `native_qfo_swiss_sequence_strata_readback_20261006.json`,
  SHA256 `8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a`.
- Original input-only Pfam inventory: `swiss_domain_annotation_inventory_20260917.json`,
  SHA256 `d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c`.
- Original bin protocol: `SWISS_DOMAIN_STRATA_PROTOCOL_20260917.md`, SHA256
  `8316866d11f988cef6a3e88f07e2802d1557a7ae581b100b85884f704b0afee8`.

Use only completed native `p0_c0_r0` and `p0_c0_r1`, with initial HMM search on.
Do not relabel these as selected high-sensitivity/satellite defaults. R0 predicts
group-clique pairs; R1 predicts resolved native pairs. R1's measurement failure
remains retained, with no inference resource admission or timing eligibility.
No reads of unfinished native outputs, parameter changes or new inference.

## Unchanged Membership And Statistics

Retain all 18 families and 563 exact accession matches, including five annotated
proteins with zero Pfam hits. Require native represented-gene membership to
agree exactly with the input-only inventory. Recheck the selected original
annotation JSONs against their recorded hashes and reproduce per-protein
features and family summaries; list unselected source hashes as inherited,
not newly checked. Coordinates remain as supplied; no width/coverage inference.

Use exactly the September 17 definitions: median distinct Pfam types below
two versus at least two (12/6 families), and repeated-type fraction below
one quarter versus at least one quarter (15/3). Higher-type families remain
APP, BAR, HOX, NOX, TRFE and VATB; higher-repeat families MAPT, PSEN and TRFE.
Include an all-family row and every family count/statistic. Refuse changed
membership or missing annotations instead of replacing bins or imputing zeros.

For each family use the benchmark's TP/FP/FN divided by two plus one prior.
Average family precision and recall within each bin, then calculate their
harmonic F1. Do not pool pairs or average family F1. Report R1-minus-R0 point
differences for F1, precision and recall, and retain both positive and negative
changes. Independently parse both native raw files and annotations, with no
import of the primary exporter, to check counts, truth, bins and exported TSV.

## Limits And Delivery

No new bootstrap draws, confidence intervals, significance tests or transfer of
historical selected-default intervals. Overlapping bins and small curated
samples do not establish independent units or causal subgroup effects. Domain
type/instance counts are not a complete validated architecture, fragment,
domain-loss or ancestral-duplication assay. Annotations share the FAS reference
resource, so do not claim independent FAS validation. No family selection or
method tuning follows these diagnostics.

Deliver machine-readable scores, family counts and a readable table plus a
separate independent readback. Record source identities, checks and execution
receipts. These are shared-host postprocessing measurements only; contention
has unknown potentially tool-dependent effects and does not establish isolated
speed rankings. Preserve original manuscript/archive/protocol bytes and all
failure evidence. Publication completion remains unproven.
