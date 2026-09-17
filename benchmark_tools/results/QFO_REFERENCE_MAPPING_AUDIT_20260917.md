# QfO Reference Mapping Audit

This audit checks the four recovered sequence-refinement stages, not every
publication competitor. It changes neither the method nor native scores.

## VGNC

The frozen reference contains 23,934 unique asserted pairs, 36,986 proteins and
16,863 family labels. Eleven proteins occur under multiple family labels.
The native scorer assigns each protein the last family label encountered in
the reference while retaining the asserted pair truth independently of that
label. Consequently, a true-positive row need not have identical family labels.
There are 24 such rows in each unrefined stage and 17 in each refined stage.

`vgnc_count_consistency_20260917.json` reconstructs native precision, recall and
harmonic F1 from TP/FP/FN counts, matching admitted endpoints within 5e-8.
All four stages have identical TP-union-FN truth sets and protein annotations.

| Stage | TP | FP | FN | Harmonic F1 |
| --- | ---: | ---: | ---: | ---: |
| multipass | 23192 | 120806 | 742 | 0.276207036 |
| multipass_refined | 19931 | 15790 | 4003 | 0.668208868 |
| strict_profiles | 23186 | 121470 | 748 | 0.275057833 |
| strict_profiles_refined | 19916 | 15806 | 4018 | 0.667694783 |

The reusable mapping audit independently reads SQLite protein mappings, maps
integer reference pairs to accession pairs, checks every raw annotation and
the exact TP/FN truth partition, and checks every emitted FP against native
cross-species family-availability rules. All four stages pass with an identical
selected-mapping digest. No TP/FP or FN/FP category overlaps occur.

The first audit rejected repeated protein numbers. Inspection of
`map_relations.py:add_reference_proteomes` showed that the native database
retains accession aliases. Allowing only identical duplicate rows was also
insufficient: protein 159158 has both Q2Z1P8 and Q2Z1P8_CANLF. Each stage's VGNC
subset has 79 extra alias rows and zero identical duplicate rows. Reconstruction
using the last rowid alias matches every retained raw label. Conflicting species
and accession assignments to different reference proteins remain fatal.
Native SQL has no explicit ordering, so this does not establish a general
alias-order guarantee. No upstream database or scorer was edited.

Reproduce from the repository root:

```bash
python benchmark_tools/audit_qfo_vgnc_mapping.py --output benchmark_tools/results/vgnc_mapping_audit_20260917.json
python -m pytest tests/unit/test_audit_qfo_vgnc_mapping.py tests/unit/test_audit_qfo_treefam_counts.py tests/unit/test_audit_qfo_swiss_counts.py -q
```

The audit hashes selected mapping content, not full prediction databases. It
does not independently rescore all predicted pairs or detect omitted FPs.
The earlier count/reference inventories retain their original, narrower scope;
the mapping report supplies subsequent evidence rather than rewriting history.
Cross-family FPs and shared proteins preclude assuming naive independent-family
resampling without additional justification. No confidence intervals were
computed and no pairwise-independent uncertainty claim is made.

## TreeFam-A Source Recovery

`treefam_source_inventory_20260917.json` records source retrievals. The QfO 2020
Zenodo record (15087752, HTTP 200) lists 25 files, including the pooled TreeFam-A
reference, but not the original NHX collection or `treefam2reference.txt`.
The successful recursive QfO repository tree query was not truncated; its
TreeFam-A JSON is an aggregation/visualization stub, not a family mapping.
The queried EBI paths returned 404 and the Sanger release path returned 403.
These observations do not prove that no recoverable archive exists elsewhere.

The inspected generator merges original trees into one native case and leaves
only the final loop tree in its serialized tree field. The completed native
count audit remains valid, but original-family uncertainty remains unresolved.
Neither pooled gene pairs nor unvalidated connected components are substituted
as independent biological families.

## Evidence

- `vgnc_reference_inventory_20260917.json`: reference identities and overlaps.
- `vgnc_raw_label_inventory_20260917.json`: raw identities and category checks.
- `vgnc_count_consistency_20260917.json`: arithmetic and cross-stage consistency.
- `vgnc_mapping_audit_20260917.json`: executable reference/database mapping audit.
- `treefam_source_inventory_20260917.json`: source recovery attempts and limits.

## Subsequent Prediction-Database Rescore

`audit_qfo_vgnc_predictions.py` now reconstructs complete native TP, FP and FN
sets from each of the four mapped prediction databases. Every category matches
the saved raw pairs exactly, not merely in count. Native P/R and harmonic F1
match admitted endpoints within 5e-8. Full databases, references and raw outputs
are checksummed before and after the audit; no input changed.

| Stage | Predictions among reference proteins | Unscored predictions |
| --- | ---: | ---: |
| multipass | 232314 | 88316 |
| multipass_refined | 40460 | 4739 |
| strict_profiles | 232521 | 87865 |
| strict_profiles_refined | 40520 | 4798 |

Unscored means neither an asserted TP nor an eligible native FP. These pairs
must not silently be added to the precision denominator as false positives.
FP classification is evaluated independently of asserted truth, matching the
native scorer; category overlap is possible in synthetic fixtures but absent
in all four retained stages. The rescore preserves native one-direction query
semantics, subset restrictions, duplicate elimination and final accession aliases.

```bash
python benchmark_tools/audit_qfo_vgnc_predictions.py --output benchmark_tools/results/vgnc_prediction_rescore_20260917.json
python -m pytest tests/unit/test_audit_qfo_vgnc_predictions.py -q
```

`vgnc_prediction_rescore_20260917.json` records exact-set digests, source identity,
full database hashes and reconstructed metrics. This supersedes the earlier
mapping-only audit's omitted-FP limitation for these four stages. It does not
audit every publication competitor, establish completeness of upstream pair
conversion, or infer independent resampling units. Remaining work includes
defensible dependent-unit uncertainty, other QfO challenges and the broader
publication requirements.
