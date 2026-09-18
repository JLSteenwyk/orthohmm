# Corrected OrthoMCL Final-Group Pair Conversion

## Semantics And Checks

`prepare_qfo_corrected_orthomcl_pairs.py` converts only independently admitted
native OrthoMCL final groups into cross-species clique pairs. It does not
use the pre-clustering similarity graph, infer additional orthologs, add
singleton groups, or silently fill ungrouped input coverage.

Before conversion it requires completed native-output admission on the
expected 2-CPU/64-GiB bizon allocation, the frozen admission executor
`fe08c7d99f266f2705e32c85294adfec676c008f`, its source identity and unchanged
admission/input/output records. The dedicated Python runtime is checked
before and after conversion. The frozen QfO assessment environment supplies
the reference mapping. The actual GG universe must still contain 984,137
proteins and 78 taxa.

The converter independently repeats the final-group/partition audit and
compares every summary field with admission. It rejects global accession
collisions, unknown or multiply grouped genes, duplicate group labels and
singleton final groups. For each group it partitions members by species,
emitting Cartesian products between distinct species exactly once. Every
pair is lexically canonical. Global output order follows native group order,
then sorted species/accessions; it is not claimed to be globally sorted.

The emitted count must agree with both the per-group species-size formula
and admission. Grouped and ungrouped protein counts must also agree. QfO
mapping must retain every emitted pair; otherwise the conversion fails and
preserves partial outputs and its failure report. Empty pair output is not
replaced with invented predictions. Downstream scoring must handle or
explicitly report such a case.

Success is `corrected_orthomcl_pairs_prepared_unscored`, with participant
`qfo_corrected_orthomcl` and semantics `cross_species_final_group_cliques`.
The original query-failure diagnostics are carried forward. Output is under
`benchmarks/results/qfo_corrected_comparator_pairs_v1/orthomcl/` and includes
`pairs.tsv`, `pairs.qfo.tsv`, `groups.json` and `results.json`.

## Verification

39 new tests cover admission status/semantics/counts/accounting, accession
collisions including ungrouped inputs, invalid or overlapping groups,
canonical cross-species expansion, missing coverage, empty predictions and
pending-job rejection. Ten deterministic random disjoint-group fixtures
match independent brute-force pair enumeration without duplicates.

The retained staged native fixture was actually tested, not skipped:
42 input proteins, 41 grouped, 12 groups and exactly30 distinct
cross-species pairs, matching brute-force enumeration. The single ungrouped
protein remains ungrouped. This is a converter fixture, not a QfO score.

Focused suite:205 passed in9.02s with legacy native runtime tests enabled.
Batch syntax passes `bash -n`. A dedicated-runtime CLI preflight against
pending admission21752 passed runtime verification and rejected the missing
completed accounting before creating an output directory. Its named
preflight allocation variables were not a real Slurm conversion run.

## Scheduling

Batch `qfo_corrected_orthomcl_pairs_batch_20260918.sh` requests2CPUs/64GiB/
24h onbizon with no requeue, clean isolated Python startup and a fresh
bytecode-cache prefix. It computes the admission manifest checksum only
after dependency satisfaction. Submit only with `afterok:21752` and a frozen
converter executor. Conversion completion still requires terminal and
artifact verification before assessment. No corrected OrthoMCL score exists
at this milestone.
