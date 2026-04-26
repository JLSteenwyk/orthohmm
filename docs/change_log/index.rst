.. _change_log:


Change log
==========

^^^^^

Major changes to OrthoHMM are summarized here.

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
