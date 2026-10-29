.. _performance:

Research and Evaluation
=======================

Manuscripts, benchmark datasets, evaluation drivers, and historical results
are maintained separately from the application checkout. The complete
pre-separation source snapshot is preserved at commit
``ebfa635ca786796287c9ff2a3caba1c2618f016e`` on the
`research-archive-20261029 branch <https://github.com/JLSteenwyk/orthohmm/tree/research-archive-20261029>`_.

* `Publication reproduction guide <https://github.com/JLSteenwyk/orthohmm/blob/ebfa635ca786796287c9ff2a3caba1c2618f016e/benchmark_tools/PUBLICATION_REPRODUCTION.md>`_
* `Manuscript and figures <https://github.com/JLSteenwyk/orthohmm/tree/ebfa635ca786796287c9ff2a3caba1c2618f016e/manuscript>`_
* `Evaluation scripts and evidence <https://github.com/JLSteenwyk/orthohmm/tree/ebfa635ca786796287c9ff2a3caba1c2618f016e/benchmark_tools>`_

OrthoBench and Quest for Orthologs are the primary benchmarks; Three
Kingdoms provides supplementary conserved-family evidence. These datasets
were inspected during development, so they do not independently validate
settings selected using their scores. Orthogroups, root hierarchical
orthogroups, and native pairwise ortholog predictions have different
semantics and must not be treated as interchangeable outputs.

Shared-host timing measurements may be distorted by competing analyses in
method-dependent ways. Matching CPU and memory limits does not establish
isolated speedups. Historical reproduction instructions and limitations
remain in the archived evidence; separating the application does not
establish new biological accuracy, timing, or reproducibility claims.
