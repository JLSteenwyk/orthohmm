Advanced Usage
==============

^^^^^

OrthoHMM output
---------------
All OrthoHMM outputs have the prefix *orthohmm* so that they are easy to find.

- orthohmm_gene_count.txt
	- A gene count matrix per taxa for each orthogroup.
|
- orthohmm_orthogroups.txt
	- Genes present in each orthogroup.
|
- orthohmm_single_copy_orthogroups.txt
	- A single-column list of single-copy orthologs.
|
- orthohmm_orthogroups
	- A directory of FASTA files wherein each file is an orthogroup.
|
- orthohmm_single_copy_orthogroups
	- A directory of FASTA files wherein each file is a single-copy ortholog.
	- Headers are modified to have taxon names come before the gene identifier.
	- Taxon names are the file name excluding the extension.
	- Taxon name and gene identifier are separated by a pipe symbol "|".
	- This aims to help streamline phylogenomic workflows wherein sequences will be concatenated downstream based on taxon names.

^^^^^

This remaining sections describe the various features and options of OrthoHMM.

- `Output directory`_
- Phmmer_
- CPU_
- `Single-copy Threshold`_
- MCL_
- `Inflation Value`_
- `All options`_

|

.. _`Output directory`:

Output directory
----------------

Output directory name to store OrthoHMM results. This directory should already exist.
By default, results files will be written to the same directory as the input
directory of FASTA files. (-o, --output_directory)

.. code-block:: shell

	# specifying output directory
	orthohmm <path_to_directory_of_FASTA_files> -o <output_directory>

.. _Phmmer:

|

Phmmer
------

Path to phmmer executable from HMMER suite. By default, phmmer
is assumed to be in the PATH variable; in other words, phmmer
can be evoked by typing `phmmer`.

.. code-block:: shell

	# specify path to phmmer executable 
	orthohmm <path_to_directory_of_FASTA_files> -p /path/to/phmmer

|

.. _CPU:

CPU
---
Number of CPU workers for multithreading during sequence search.
This argument is used by phmmer during all-vs-all comparisons.
By default, the number of CPUs available will be auto-detected.

.. code-block:: shell

	# run orthohmm using 8 CPUs 
	orthohmm <path_to_directory_of_FASTA_files> -c 8

|

.. _`Single-copy Threshold`:

Single-copy Threshold
---------------------

Taxon occupancy threshold when identifying single-copy orthologs.
By default, the threshold is 50% taxon occupancy, which is specified
as a fraction - that is, 0.5.

.. code-block:: shell

	# specify single-copy threshold as a fraction 
	orthohmm <path_to_directory_of_FASTA_files> -s 0.5

|

.. _MCL:

MCL
---

Path to mcl executable from MCL software. By default, mcl
is assumed to be in the PATH variable; in other words,
mcl can be evoked by typing `mcl`.

.. code-block:: shell

	# specify path to mcl executable 
	orthohmm <path_to_directory_of_FASTA_files> -m /path/to/mcl

|


.. _`Inflation Value`:

Inflation Value
---------------

MCL inflation parameter for clustering genes into orthologous groups.
Lower values are more permissive resulting in larger orthogroups.
Higher values are stricter resulting in smaller orthogroups.
The default value is 1.5.

.. code-block:: shell

	# use an inflation value of 1.5 during mcl clustering 
	orthohmm <path_to_directory_of_FASTA_files> -i 1.5

|

.. _`All options`:

All options
---------------------


+------------------------------+--------------------------------------------------------------------------------+
| Option                       | Usage and meaning                                                              |
+==============================+================================================================================+
| -h/\-\-help                  | Print help message                                                             |
+------------------------------+--------------------------------------------------------------------------------+
| -v/\-\-version               | Print software version                                                         |
+------------------------------+--------------------------------------------------------------------------------+
| -o/\-\-output_directory      | Output directory name. Default: same directory as directory of FASTA files     |
+------------------------------+--------------------------------------------------------------------------------+
| -p/\-\-phhmer                | Path to phmmer from HMMER suite. Default: phmmer                               |
+------------------------------+--------------------------------------------------------------------------------+
| -c\-\-cpu                    | Number of parallel CPU workers to use for multithreading. Default: auto detect |
+------------------------------+--------------------------------------------------------------------------------+
| -s/\-\-single_copy_threshold | Taxon occupancy threshold for single-copy orthologs. Default 0.5               |
+------------------------------+--------------------------------------------------------------------------------+
| -m/\-\-mcl                   | Path to mcl software. Default: mcl                                             |
+------------------------------+--------------------------------------------------------------------------------+
| -i/\-\-inflation_value       | MCL inflation parameter. Default: 1.5                                          |
+------------------------------+--------------------------------------------------------------------------------+

