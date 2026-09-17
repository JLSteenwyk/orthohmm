.. _performance:


Performance Assessment
======================

OrthoHMM is evaluated against other orthology-inference tools using OrthoBench
and Quest for Orthologs (QfO) as primary benchmarks. Three Kingdoms is a
supplementary conserved-family benchmark. Results depend on the dataset,
configuration, reference universe and type of prediction being evaluated.
The current publication analyses do not establish general superiority over
full OrthoFinder or the other comparators.

Current evidence
----------------

* On OrthoBench, phylogenetic OrthoHMM has a higher observed aggregate F1 than
  full OrthoFinder, but the paired difference interval includes zero. The
  results show a precision-recall trade-off, not established F1 superiority.
* QfO challenges measure different properties. Their individual endpoints
  must remain visible; the project's six-metric mean is a secondary summary,
  not an official overall accuracy score.
* Three Kingdoms scoring is restricted to the available BUSCO reference
  universe. It does not establish proteome-wide orthology accuracy.
* Frozen YGOB evaluation provides bounded transfer evidence on novel taxa,
  but substantial family overlap prevents a family-disjoint validation claim.
* The prespecified biological application shows improved homolog-supported
  paralog separation with phylogenetic OrthoHMM relative to high sensitivity,
  alongside lower homolog coverage. Full OrthoFinder and SonicParanoid perform
  better on supported separation and coverage in this application.
* Matched-resource timing runs are underway on a dedicated machine. Historical
  resource measurements should not be pooled with this new timing panel.

OrthoBench and QfO were inspected during development. Their scores are not
independent confirmation of settings selected using those datasets. The
repository retains neutral and negative results, native failures, exclusions,
paired uncertainty and output-conversion audits.

Prediction semantics
--------------------

Orthogroups, root hierarchical orthogroups and native pairwise orthologs are
different outputs. Pre-clustering graph edges are not final ortholog predictions.
The retained OrthoFinder MCL checkpoint is a diagnostic comparison, not a
separately executed sequence-only OrthoFinder pipeline. Comparisons must use
the output and scoring protocol specified for each benchmark.

Reports and reproducibility
---------------------------

The `eight-method comparison <https://github.com/JLSteenwyk/orthohmm/blob/main/benchmark_tools/results/PUBLICATION_COMPARISON_ORTHOMCL_COMPLETE_20260916.md>`_
contains the consolidated observed scores. The
`claim-to-evidence checklist <https://github.com/JLSteenwyk/orthohmm/blob/main/benchmark_tools/results/PUBLICATION_CLAIMS_20260916.md>`_
links uncertainty, ablations, simulations, transfer evaluation and the
biological application. The
`progress ledger <https://github.com/JLSteenwyk/orthohmm/blob/main/benchmark_tools/results/PUBLICATION_PROGRESS.md>`_
records completed analyses, running jobs and unresolved requirements.

The new publication package remains in preparation; its release and archival
deposition are not complete. Earlier desirability rankings and the historical
performance illustration remain in the repository and its history, and are
not substituted for the current audited endpoints.
