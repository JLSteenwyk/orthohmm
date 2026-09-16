# OrthoMCL BLAST Failure Impact

This audit assesses direct reference exposure of the 53 failed queries,
not the counterfactual performance of a repaired search and clustering run.
The prior input/log/group audit establishes that all 53 are absent from
final OrthoMCL groups. They represent 0.0054275% of the 976,504 input proteins.

## Direct Exposure

| QfO endpoint | Failed-query exposure | Interpretation |
| --- | --- | --- |
| SwissTrees | No failed protein in its 18 mapped reference cases | No reference relation directly incident to these failed proteins |
| TreeFam-A | No failed protein in its mapped reference case | No reference relation directly incident to these failed proteins |
| VGNC | Zero incident pairs among 23,934 asserted pairs | No directly missing asserted pair attributable to these proteins |
| EC | Zero failed proteins with EC annotations | No direct loss of an EC-annotated pair endpoint |
| GO, experimental evidence | Four failed proteins with qualifying annotations | Potential direct coverage loss; not a measured GO-score change |
| FAS | All 53 have entries; 46 have at least one feature type | Feature presence alone does not establish pair eligibility or score impact |

The four experimental-GO proteins are P0DKJ0 (5 amino acids, three terms),
P62120 (25 aa, one term), P62945 (25 aa, seven terms), and Q962S2 (25 aa,
three terms). The evidence filter is EXP/IDA/IPI/IMP/IGI/IEP, matching the
reviewed QfO configuration. Per-protein FAS feature-type counts are retained
in the machine-readable result; empty records must not be described as
feature-bearing proteins.

TreeFam-A is encoded as one pooled native `RecTreeCase` with 11,140 mapped
proteins. The case count is not a claim that it contains only one biological
gene family. The audit uses the benchmark's native Darwin structures, not
accession substring matching or a custom parser of its reference program.

## Baseline Decision

Retain the unmodified OrthoMCL 1.4 run as the primary comparator and disclose
the failed queries. Do not silently disable low-complexity filtering, edit
short sequences, replace the search engine, remove genes from the benchmark
universe, or impute orthology. Such changes would create a different run or
method configuration, not repair the provenance of this baseline.

The observed reference exposure does not justify another full BLAST run to
replace this baseline. A selective unmasked-search sensitivity experiment
could test the statistics failures, but must retain the original result and
be labeled a modified-search diagnostic. No such experiment has been run,
and the short-query failures may remain unsupported by the legacy engine.

No zero-effect or negligible-effect claim follows from these counts. New
edges can change clusters containing other genes; direct reference absence
does not bound those indirect effects. Functional benchmarks score selected
predicted pairs rather than a complete truth set. Neither a change in their
mean similarity nor a precision/recall repair benefit has been measured.

## Reproduction And Evidence

`audit_orthomcl_reference_impact.py` invokes the same Darwin container used by
QfO, reads mapped proteins and reference relations, checks VGNC unordered
pairs, and records input, source, container, annotation, and native-log hashes.
`orthomcl_reference_impact.drw` mirrors the native experimental-GO evidence
filter and reference-case size threshold. The result is
`orthomcl_reference_impact_20260916.json`.

Nine new tests cover native-report completion, mapped identifiers, relation
orientation, error detection, annotation counts, and unique VGNC pairs.
Together with the original BLAST audit tests, 17 focused tests pass. Native
execution also completed without reported errors or warnings. Raw audit
outputs remain in `benchmarks/work/orthomcl_reference_impact_v2`; the first
audit is preserved in `orthomcl_reference_impact_v1` and was superseded to
clarify native case counts and distinguish FAS entries from feature content.
