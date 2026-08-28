.. _change_log:


Change log
==========

^^^^^

Major changes to OrthoHMM are summarized here.

**0.5.0**
Sets Leiden CPM at resolution 0.1 as the default clustering path after
head-to-head benchmarking across OrthoBench, QfO 2020, and Three Kingdoms.
The default fixed setting is best for the highest-priority OrthoBench and
QfO benchmarks. ``--cpm_resolution auto`` remains available for distant
cross-kingdom inputs, where it improves Three Kingdoms substantially
(F=0.907 vs 0.0003 for fixed 0.1). ``--clustering mcl`` remains available
as an optional conservative fallback, but mcl is no longer required for
the default pipeline.

**0.4.2**
``--cpm_resolution auto`` now anchors γ on the smallest *positive* edge
weight (excluding numerical-zero artifacts). The earlier formula
(γ = 4 × strict-min) collapsed to γ ≈ 0 on the QfO 40M-edge graph and
caused Leiden to segfault.

**0.4.1**
Refines ``--cpm_resolution auto`` to use γ = 4 × min(edge_weight),
beating the v0.4.0 10th-percentile heuristic on both OrthoBench (F=65.9
vs 58.8) and Three Kingdoms (F=0.907 vs 0.821).

**0.4.0**
New ``--cpm_resolution auto`` flag. Auto-tunes γ to the post-RBNH
edge-weight distribution so OrthoHMM no longer collapses to
all-singletons on distantly-related inputs.

**0.3.0**
Replaced the external MCL binary with in-process Leiden CPM clustering
via ``igraph`` + ``leidenalg``. Both libraries are pip-installed wheels,
so OrthoHMM no longer requires any external executable when run with
defaults. Leiden CPM (resolution=0.1) beats MCL on the OrthoBench 2020
reference: F=65.7% vs MCL's best F=62.4% (inflation=1.5) on the
identical RBNH edge set. New flags: ``--clustering {leiden, mcl}`` and
``--cpm_resolution``. Selecting ``--clustering mcl`` reverts to the
prior MCL pipeline and re-introduces the external ``mcl`` requirement.

**0.2.0**
Added a built-in profile HMM + k-mer prefilter search engine that
replaces the ``phmmer`` subprocess. The new engine is the default and
substantially reduces wall time and memory on multi-proteome datasets
(see the bacterial scaling table on the home page). HMMER is now
optional — only required when opting into ``--search_mode phmmer``.
Also adds the WAG and LG substitution matrices, drops Python 3.9
support (now requires Python 3.10+), and adds optional C/AVX2 + CUDA
kernels that are compiled at install time when a suitable toolchain
is available; otherwise the runtime falls back to a Numba
implementation transparently.

**0.1.1**
There is no longer a limit on the length of gene names for single-copy orthologous genes.

**0.1.0**
Modified how to handle phmmer multiprocessing, giving the user a parallelized experience.
Specifically, if a user sets CPUs to 8, 8 runs of phmmer will run at the same time.
