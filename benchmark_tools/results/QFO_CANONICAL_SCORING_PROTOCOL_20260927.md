# Canonical QfO Scoring Protocol

## Scope Frozen Before Scoring

Evaluate the validated canonical-order predictions from native 22333 and
independent readback 22334. Do not tune parameters, choose between retained
and canonical order based on scores, change endpoint definitions, or overwrite
historical results. This is a development-exposed ordering/reproducibility
analysis, not independent generalization or a new method-selection experiment.

Native plan SHA256:
`465a0509d7d00640d32cde8411b145ad454cf52627be8469083e37f695fb195d`.
Validated three-way result SHA256:
`e1409d69b48565ff482a70be83b0e7859ec03be896ebbc015489bd699f6c0ab8`.
Require successful scheduler states and the bound final readback execution
receipt; recheck native admission and recorded artifact identities before
conversion. No additional inference or tree rerun is needed.

## Conversion

Use the existing `prepare_qfo_factorial_pairs.write_native_pairs` converter,
`simulation_method_outputs.orthohmm_pairs` parser, accession normalization and
`qfo_filter_pairs` reference mapping. Input is the native pairwise ortholog TSV,
not root-group cliques, graph edges or confidence-weighted pairs. Check sequence
ownership against the admitted 78 FASTAs/provenance and require unique accession
normalization, unique ordered native pairs and cross-species predictions.

The expected native total is 5,959,535 pairs. Require this exact converted count
and zero reference-mapping loss, as in the corrected historical comparison.
Unexpected loss or duplicates must stop admission, not silently modify the
evaluated universe. Preserve raw conversion, mapped output, counts, checksums,
source identities, scheduler record and logs in a fresh canonical namespace.
Use a scheduled 2-CPU conversion; do not report it as inference runtime.

## Frozen Scoring Environment

Use `qfo_assessment_environment_20260917.json`, SHA256
`e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc`.
All 704 records returned by `run_qfo_recovered_assessment.environment_records`
were rehashed successfully before freezing this protocol: pipeline sources,
references, containers, Java, executables and support/configuration files.
The existing `command_for` helper provides the native Nextflow invocation.

- Participant: `ohmm_qfo_canonical_20260927`.
- Event year: 2020.
- Challenges: GO, EC, VGNC, SwissTrees, TreeFam-A, FAS; retain all six.
- References: the frozen pipeline's `reference_data/2020` and
  `reference_data/data`, including the same `mapping.json.gz` as conversion.
- Work directory: `qfo_benchmark/w/qcan27` (within the validated Darwin path limit).
- Native score outputs: `qfo_benchmark/scoring/canonical_20260927`.
- Receipts/logs: `benchmarks/results/qfo_canonical_assessment_20260927`.
- Conversion outputs: `benchmarks/results/qfo_canonical_pairs_20260927`.
- Scoring resources: 8 CPUs, 128 GiB RAM, 12 hours, one attempt, no requeue.

Keep the original execution configuration, offline Nextflow, pinned Java and
container cache settings. No `-resume` into a prior scoring namespace. Freeze
the exact command, converted-input identity and executable source records before
submission. A successful process alone is insufficient: independently validate
native assessment artifacts, challenge coverage, participant identity and
prediction counts using existing assessment validators before exporting scores.
Retain failure diagnostics and do not automatically retry a failed scorer.

## Reporting And Interpretation

Report each endpoint, its native coverage/prediction count, available precision
and recall, and canonical-minus-historical differences under unchanged scoring
definitions. Use the corrected historical `p1_c1_r1` arm as comparator only after
rechecking its admitted score provenance. Preserve historical scores as a separate
row. The six-endpoint mean remains a project-defined secondary summary, not F1
or a universal ranking. No new independent confidence interval is implied.

The frozen FAS scorer calls `random.shuffle` without a fixed seed. Preserve
the new native sample and its artifacts; label FAS differences from the saved
historical sample as confounded by sampling variability. Do not attribute such
a difference solely to ordering, relabel the native SEM as a paired interval,
or repeatedly rescore until a preferred value appears. Any later matched-sample
analysis needs its own prespecified protocol and must retain native results.

Keep all negative/neutral outcomes. The already established 110 removed pairs,
85 added pairs and seven annotation changes are prediction differences, not
accuracy changes. This evaluation does not close remaining generalization,
uncertainty, dedicated timing, biological-validation or publication-release gaps.
