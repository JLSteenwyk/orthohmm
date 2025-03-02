.. _change_log:


Change log
==========

^^^^^

Major changes to OrthoHMM are summarized here.

**0.1.1**
There is no longer a limit on the length of gene names for single-copy orthologous genes.

**0.1.0**
Modified how to handle phmmer multiprocessing, giving the user a parallelized experience.
Specifically, if a user sets CPUs to 8, 8 runs of phmmer will run at the same time.
