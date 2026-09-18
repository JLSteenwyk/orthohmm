# Corrected Proteinortho Assessment

Corrected-input Proteinortho 6.3.6 assessment **21718** completed successfully
in scheduler elapsed **36:06** (8 CPUs, bizon). Independent admission
**21719** completed in **17 seconds** (2 CPUs). These are assessment and
validation durations, not inference runtime or dedicated scaling measurements.

The unchanged frozen validator was run again to
`benchmarks/work/qfo_corrected_proteinortho_assessment_readmission_20260918.json`.
Its entire parsed report agrees exactly with the first admission, including
all 48 native assessment records, 15 successful native tasks and 14 metric files.
The retained report is `qfo_corrected_proteinortho_assessment_20260918.json`,
SHA-256 `00acc10285b955cf0b33c8b3d95044e4f797ab77ba6ae4e768b62615e61f1552`.

## Endpoints

Values below are the report's native endpoints, not six interchangeable F1s.

| Endpoint | Score | Recall | Precision |
| --- | ---: | ---: | ---: |
| GO average Schlicker | 0.486336400 | n/a | n/a |
| EC average Schlicker | 0.963168110 | n/a | n/a |
| VGNC F1 | 0.954896244 | 0.930559037 | 0.980540636 |
| SwissTrees F1 | 0.718111000 | 0.576414360 | 0.952179940 |
| TreeFam-A F1 | 0.643187434 | 0.487117070 | 0.946415030 |
| FAS | 0.813594864 | n/a | n/a |

The project-defined secondary six-metric mean is **0.7632156753824589**.
It is not an official QfO F1 or a basis for claiming general superiority.

The submission contains **4,695,385** native post-clustering cross-species
pairs; all mapped to the corrected reference. GO assesses 89,609 relations,
EC 143,135; these challenge counts are not the total submitted coverage.
Input universe: 78 proteomes, 984,137 unique sequences. This result must not
replace the original-release row or be compared with old-release scores as
though the input data were identical.

## Provenance And Remaining Work

- Scoring executor: `74afad5376b7ee11fdabfba386851fd8d3c02857`.
- Independent validator: `f7e80d3a94cc805f7a09c646c50a6c5c4a656343`.
- Conversion job: 21717; pair manifest SHA-256
  `c0620eb8772a983956d2f3b5cb98636bbaca2ba0c2adbcd9ee0ff0805d02342b`.
- Participant: `qfo_corrected_proteinortho`; frozen 2020 QfO assessment
  environment unchanged. See the native/conversion and assessment-submission
  ledgers for inference inputs, graph semantics and exact commands.

Paired corrected-release family uncertainty has not been computed. Native
error bars are not paired method-difference intervals; GO/EC's fields labeled
`stderr` contain the previously audited native confidence half-widths.
TreeFam family-source uncertainty remains unresolved. FAS includes its native
unseeded sampling, so a new assessment could differ even with identical inputs.
Other corrected comparator rows and the factorial experiment remain pending;
no corrected-release ranking or publication-readiness claim is supported here.
