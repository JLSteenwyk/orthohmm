.. image:: _static/img/logo_condensed.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/orthohmm

OrthoHMM
========

OrthoHMM infers orthogroups using a built-in profile-HMM search engine and
Leiden clustering. The standard pipeline does not require external
bioinformatics executables. Optional modes support high-sensitivity profile
expansion, phylogenetic reconciliation, and legacy HMMER/MCL inference.

Installation
------------

Use Python 3.10 or newer and an isolated environment. These commands install
the development application from source, not a potentially older PyPI release.

.. code-block:: shell

   git clone --depth 1 https://github.com/JLSteenwyk/orthohmm.git
   cd orthohmm
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install .
   orthohmm --help

A C compiler with OpenMP support builds the native CPU kernels during
installation. On Linux, install GCC and its development tools. On macOS,
install Homebrew GCC and ensure ``gcc`` resolves to that compiler when
installing. Standard search can fall back to Numba when compilation is
unavailable; high-sensitivity MSA profile expansion requires native kernels.
Set ``ORTHOHMM_CPU_TARGET=baseline`` during installation when building for
compatible machines other than the build host.

Usage
-----

Place one protein FASTA file per species in an input directory and create the
output directories before running.

.. code-block:: shell

   mkdir results sensitive_results
   orthohmm proteomes/ -o results/ -c 8
   orthohmm proteomes/ -o sensitive_results/ -c 8 --accuracy_profile high_sensitivity

See :doc:`advanced/index` for options and output files.

Optional dependencies
---------------------

* Phylogenetic reconciliation: install ``'.[phylogeny]'`` and provide MAFFT
  and FastTree executables. A supplied or inferred species tree can be used.
* Legacy sequence search: install HMMER and select ``--search_mode phmmer``.
* Legacy clustering: install MCL and select ``--clustering mcl``.
* Python HMMER control helpers: install ``'.[hmmer]'`` for their Biopython
  dependency. This extra is not needed by the built-in search engine.

Research and citation
---------------------

The application repository keeps source, tests, and user documentation.
Manuscripts, benchmark datasets, campaign scripts, and historical results
are preserved on the
`research-archive-20261029 branch <https://github.com/JLSteenwyk/orthohmm/tree/research-archive-20261029>`_.
See :ref:`performance` for the research snapshot and evaluation limitations.

Please cite *OrthoHMM: Improved Inference of Ortholog Groups using Hidden
Markov Models*, Steenwyk et al. (2024), bioRxiv,
`doi:10.1101/2024.12.07.627370 <https://doi.org/10.1101/2024.12.07.627370>`_.
Citation metadata are also available in ``CITATION.cff``.

.. toctree::
   :maxdepth: 4

   about/index
   advanced/index
   performance_assessment/index
   change_log/index
   other_software/index
   frequently_asked_questions/index
