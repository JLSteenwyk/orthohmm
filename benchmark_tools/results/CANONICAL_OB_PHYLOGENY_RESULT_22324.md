# Fresh Canonical OrthoBench Phylogeny Result

Native job 22324 completed 0:0 in 25:58 with 32 allocated CPUs. Its pinned
post-run audit, job 22325, completed 0:0 in 2:30 with two CPUs. Both used
single attempts and no requeue. The
[result receipt](canonical_ob_phylogeny_result_22324.json) binds scheduler
accounting, frozen plans, submissions, native completion, all four scientific
readbacks and full score comparisons. These records and the 714 pinned local
Python sources were rechecked after terminal completion.

## Partition and Accuracy

| Run | Root groups | Precision (%) | Recall (%) | F1 (%) |
| --- | ---: | ---: | ---: | ---: |
| Historical baseline | 59,770 | 81.770454 | 67.755336 | 74.106074 |
| Fresh installed full run | 59,812 | 81.103625 | 67.739444 | 73.821569 |
| Fresh canonical-candidate phylogeny | 59,770 | 81.770454 | 67.755336 | 74.106074 |

The new root partition equals the historical partition exactly after ignoring
group labels: all 59,770 groups match, with no changed gene memberships.
The complete frozen score objects also match, including all 70 reference-family
records and 15 exactly recovered RefOGs. No historical score is replaced.

Relative to the separately installed full run, F1 differs by +0.2845044525508911
percentage points. There are 26,433 genes in changed groups and eight changed
RefOG records (007, 011, 021, 023, 027, 035, 053, 058), all retained in the
receipt. This is a reproducibility comparison, not a newly independent accuracy
test or an adjusted significance claim.

## Scope and Validation

The experiment used validated fresh search evidence, the private Leiden 0.11
distribution and experimental canonical directed-pair ordering to recover
historical candidate inputs. It then inferred all alignments and species/gene
trees afresh with the frozen rules and membership constraints. It did not
reuse old phylogeny checkpoints or repeat the expensive all-to-all search.

Independent readers checked structure, sequence content, reconciliation
events/pairs and hierarchy semantics. The run covers all 251,378 proteins in
12 species, with 54,445 candidate families, 8,681 reconciled families, 45,764
bypassed families, 200 species-tree families and 966,439 native ortholog pairs.
All gene-tree and species-tree checkpoint reuse counters are zero.

This establishes recovery of the historical root partition and score in this
bounded OrthoBench experiment. It does not prove general determinism across
machines or dependencies, biological tree truth, QfO reproduction, or separate
causal attribution of final F1 to the dependency versus ordering change.
The experimental ordering policy is still not a production default, and the
publication installation/reproduction workflow requires final integration.

## Resources and Remaining Work

The replay records 1,555.403741 wall seconds, 43,005.304147 user CPU seconds,
3,227.046545 system CPU seconds and 1,693,839,360 bytes sampled peak summed
process-tree RSS. These are shared-host, phylogeny-stage measurements, not full
pipeline or controlled comparative timing. The malformed live Slurm AveCPU
observation was excluded and is not used to repair or replace measured CPU time.

The protocol and tests were committed before evaluation. Prior focused suites
passed 115 workflow/reader tests and 98 launcher/pipeline/constraint tests;
the new evidence is the completed native run and independent full-data audits.
No scientific implementation was changed during inference or readback.
The full publication goal remains open, including generalization limitations,
other uncertainty endpoints, controlled scaling, remaining provenance and
release/archive packaging. DGX remains deferred.
