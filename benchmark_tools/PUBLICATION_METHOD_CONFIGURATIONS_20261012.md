# Retained Method Configurations

This reader catalogue consolidates recorded settings and command scopes. It
does not create historical argv, alter predictions, or launch a new benchmark.
Exact paths, argv arrays, input/output hashes and measurement records remain in
the [24-cell register](results/all_benchmark_metadata_integrated_20261004_v3/register.json).
Read nested metadata and supplemental fields as well as the base `commands`
field: a null or empty base field does not mean every launch detail is absent.

## Scientific Source And Output

The prospective publication method is OrthoHMM 0.5.0 at
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`, not current development HEAD.
The [YGOB protocol](results/YGOB_VALIDATION_PROTOCOL_20260916.md) freezes
high sensitivity with Leiden seed 4, BLOSUM62, E-value 1e-4, CPM resolution
0.1, default refinement, 32 CPUs and eight threads per worker. It sets
`PYTHONHASHSEED=0` and limits BLAS/OpenMP environment threads to one.
Satellite_v2 adds expanded candidates, an internally inferred species tree,
min_variance species-tree rooting, species_overlap root rule and
positive_paralogy pair rule. Its YGOB endpoint uses final root HOGs for
homolog-group recovery, not resolved pairwise orthology.

Historical rows declare versions; they do not all establish execution from
this later scientific source. The [frozen scoring specification](results/FROZEN_HMM_SCORING_SPECIFICATION_20260930.md)
and [source/runtime guide](PUBLICATION_SOURCE_COMPONENT.md) distinguish source,
workflow and native build overlays. Newer native arms retain their own
[terminal admissions/failures](results/native_qfo_terminal_failures_20261009_v1/report.json).
Do not substitute a historical cached result, a later native arm or current
HEAD for a pinned experiment merely because its method label matches.

Across the register the retained version declarations are OrthoHMM 0.5.0,
OrthoFinder 3.1.5, SonicParanoid 2.0.9, ProteinOrtho 6.3.6, FastOMA 0.3.5 and
OrthoMCL 1.4. Eight method rows include two OrthoHMM modes and two OrthoFinder
outputs. Output semantics are benchmark-specific: native pairs, group-derived
cross-species pairs, final groups and diagnostic checkpoints are not interchangeable.
The [scoring audit](results/SCORING_AUDIT_20260910.md) records conversion corrections.

## OrthoBench

Settings below are retained observations, not instructions to repeat completed
inference. The register keys match the first column.

| Method Key | Recorded Configuration And Execution Scope |
| --- | --- |
| `orthohmm_high_sensitivity` | Checked cached-hit replay: `replay_high_sensitivity.py`, `hits_BLOSUM62_mc100.pkl`, `--cpu 32 --leiden-seed 4`. This argv is replay plus official scoring, not the original initial HMM search or its full cost |
| `orthohmm_phylogeny_satellite_v2` | Native entrypoint records `-c 32 --threads_per_worker 8 -x BLOSUM62 -e 0.0001 --clustering leiden --cpm_resolution 0.1 --refinement_profile default --accuracy_profile high_sensitivity --stop infer`; reconciliation uses `satellite_v2`, inferred species tree, `species_overlap`, `positive_paralogy`, `min_variance`, and retained MAFFT 7.525/FastTree v220 paths |
| `orthofinder_3_1_5_full` | Native `orthofinder -f INPUT -t 32 -a 8 -S diamond`; omitted flags are not newly specified here |
| `orthofinder_3_1_5_sequence_only` | `orthofinder_mcl_to_orthogroups.py` converts the same full run's `clusters_OrthoFinder_I1.2.txt_id_pairs.txt` with `SequenceIDs.txt`. Not a separate inference run or independently measured sequence-only pipeline |
| `sonicparanoid_2_0_9` | Original argv unavailable. Retained native settings report default mode, DIAMOND very-sensitive, 32 threads, minimum bitscore 40, in-paralog length difference 0.75, table-merging threshold 0.75 and MCL inflation 1.5. Tool-reported duration exists; CPU time and peak RSS do not |
| `proteinortho_6_3_6` | Base argv is null, but nested native metadata retains two `proteinortho` calls on the 12 named FASTAs with `-cpus=32 -project=orthobench`. Not proof of one uninterrupted inference interval. Input sanitization changed 177 proteins/deleted 869 asterisks; its accuracy effect was not measured |
| `fastoma_0_3_5` | Base argv is null, but the retained Nextflow log records `nextflow run FastOMA.nf -profile docker`, input/output folders and work directory. Supplied species tree; final orthologous groups scored, not root-HOG diagnostic. Ten exit-137 tasks remain part of historical workflow accounting; image/database identities and peak RSS unresolved |
| `orthomcl_1_4` | Native mode 1, BLAST 2.2.13 `blastp -e 1e-5 -m 8 -a 8 -v 1000 -b 1000`, `formatdb -p t`, and MCL 02-063 `-I 1.5`. Preserve those historical eight-thread BLAST settings, not a reconstructed 32-thread command |

The SonicParanoid settings are in
`retained_ob_metadata.native_metadata.settings`. A bounded launcher-history
search found no matching source changes for its April 20 run in the selected
January-April 2026 paths. That is not global proof of absence; original argv
stays null. ProteinOrtho calls are in `native_metadata.system_calls` and
FastOMA's log lines in `native_metadata`, under their respective retained
OrthoBench metadata. No private shell/session history was consulted.

## Corrected QfO

The eight-method corrected comparator snapshot is separate from later native
component experiments. Its [score/admission register](results/qfo_corrected_comparison_20260926_v7/manifest.json)
and [selected OrthoHMM stage addendum](results/qfo_orthohmm_stage_metadata_20261004/register.json)
provide the retained identity links; empty legacy command arrays are preserved.

| Method Key | Recorded Configuration And Execution Scope |
| --- | --- |
| `orthohmm_high_sensitivity` | Selected cached cell `p1_c0_r0`: checked strict-profile seed, no candidate expansion/reconciliation, group-derived cross-species pairs. Stage addendum supplies replay/preparation metadata, not original initial-search argv or an exact selected full-pipeline cost |
| `orthohmm_phylogeny_satellite_v2` | Selected cached cell `p1_c1_r1`: same strict seed, expanded candidates and native reconciliation pairs. Replay is shared with the HS row, not an independent repetition. Retained preparation/reconciliation stages do not supply a full-pipeline cost |
| `orthofinder_3_1_5_full` | Native `orthofinder -f INPUT -t 32 -a 32 -S diamond`; full phylogenetic output |
| `orthofinder_3_1_5_sequence_only` | Register records the shared full-run command above. Prediction is its pre-phylogenetic MCL checkpoint, not another invocation or independent full-run time |
| `sonicparanoid_2_0_9` | Native `sonicparanoid -i INPUT -o OUTPUT -t 32` |
| `proteinortho_6_3_6` | Singularity `exec --bind /mnt` with retained 6.3.6 image; `proteinortho -project=qfo -cpus=32` and the complete 78-file list in the register |
| `fastoma_0_3_5` | Nextflow `-C fastoma_corrected_execution.config`, `FastOMA.nf`, `LUCA.h5`, `--fasta_header_id_transformer UniProt --force_pairwise_ortholog_generation true --max_cpus 180 --max_memory 700.GB`. These are workflow ceilings, not measured usage or matched 32-CPU resource evidence |
| `orthomcl_1_4` | Recovered native Perl mode 4 using retained `all.bpo`/`all.gg` through the wrapper. Recorded downstream interval excludes BLAST and BLAST-to-BPO conversion. All 53 sequence-specific failed queries and the [reference-impact analysis](results/ORTHOMCL_FAILURE_IMPACT_20260916.md) remain explicit |

The [stage interpretation](results/QFO_ORTHOHMM_STAGE_LINKAGE_RESULT_20261004.md)
labels shared cached replay, candidate preparation and reconciliation separately.
Do not sum their walls or peaks, impute initial search, or use these historical
scopes as isolated speed comparisons.

## Supplementary Three Kingdoms

Seven historical rows have empty base `commands`; their
`supplemental_commands` retain metrics/harness/GNU-time evidence with distinct
scope labels. The matched SonicParanoid row has its native command in the base
field. The source register preserves every complete argv and evidence pin.

| Method Key | Recorded Configuration And Execution Scope |
| --- | --- |
| `orthohmm_high_sensitivity` | Metrics entrypoint and harness launch record the 32-CPU/eight-thread, BLOSUM62/1e-4/Leiden/CPM 0.1/default-refinement/high-sensitivity settings, `--stop infer`; no phylogeny. Harness manifest was collected after execution |
| `orthohmm_phylogeny_satellite_v2` | Same base settings plus reconciliation, inferred species tree, satellite_v2, species_overlap, positive_paralogy and min_variance; retained MAFFT/FastTree paths. Metrics entrypoint and after-execution harness launch are separately labeled |
| `orthofinder_3_1_5_full` | Historical native command `orthofinder -f INPUT -t 32 -a 32 -S diamond` |
| `orthofinder_3_1_5_sequence_only` | Same-run MCL/SequenceIDs conversion only; not separate sequence-only inference |
| `sonicparanoid_2_0_9` | Matched-input native `sonicparanoid -i INPUT -o OUTPUT -t 32`, not the replaced older unmatched row |
| `proteinortho_6_3_6` | Singularity 6.3.6 image, `proteinortho -project=tk -cpus=32` and the 12 named FASTAs |
| `fastoma_0_3_5` | Nextflow Docker profile, `--container_version 0.3.5`, `LUCA.h5`, `noop` FASTA transformer, forced pairwise generation, `--max_cpus 128 --max_memory 700.GB`. Supplied-tree workflow; GNU-time CPU/RSS exclude aggregate Docker tasks. Ceilings are not actual allocation/peak usage |
| `orthomcl_1_4` | Retained BLAST-to-BPO converter followed by native mode 4 on `all.accelerated.bpo`/`all.gg`, `ORTHOMCL_PAIR_WORKERS=11`. Separate conversion/downstream intervals do not establish full BLAST-plus-inference costs |

## Resources And Reproduction

Only the separately frozen
[27-attempt Threadripper panel](results/threadripper_shared_panel_snapshot_20261004_v27/panel.json)
provides the designated matched-limit scaling experiment. It retains 25
measured/24 eligible observations and failed/ineligible identities 0, 17, 20;
there are no automatic retries. Historical commands above are not that panel.
Shared-host contention has unknown, potentially tool-dependent effects, so
matching limits does not establish isolated performance or a speed ranking.

Use the [current study delivery index](PUBLICATION_STUDY_DELIVERY_20261012.md)
for copied executable main-result commands, runtime locks and the versioned
local archive. Unknown historical consumption, argv or costs remain visible
limitations; this catalogue is consolidation, not all-tool hermetic certification.
