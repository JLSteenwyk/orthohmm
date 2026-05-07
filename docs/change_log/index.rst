.. _change_log:


Change log
==========

^^^^^

Major changes to OrthoHMM are summarized here.

**0.5.0**
Reverted the v0.3.0 default clustering choice. After head-to-head
benchmarking across phylogenetic-diversity ranges (OrthoBench 12
bilaterians, QfO 2020 78 mixed-domain proteomes, Three Kingdoms
12 cross-kingdom proteomes), MCL with inflation=1.5 is the only
clustering choice that never catastrophically fails: Leiden CPM with
the fixed cpm=0.1 default collapses to all-singletons on cross-kingdom
inputs (F=0.0003 on Three Kingdoms), and Leiden with --cpm_resolution
auto over-merges on QfO functional metrics. MCL@1.5 trails the best
Leiden setting by ~3 pp on closely-related inputs but is robust across
the full diversity range. ``--clustering leiden --cpm_resolution auto``
remains the recommended choice for users who prefer the pure-Python
path or who can confirm their inputs are tightly clustered. mcl is
once again a required external binary; see the README/install docs.

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
