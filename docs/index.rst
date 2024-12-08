.. image:: _static/img/logo_condensed.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/orthohmm

^^^^^


OrthoHMM using high sensitivity and specificity Hidden Markov Models for orthology inference.

If you found OrthoHMM useful, please cite *MANUSCRIPT TITLE*. Steenwyk et al. 2020, JOURNAL. doi: |doiLink|_.

.. _doiLink: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007
.. |doiLink| replace:: 10.1371/journal.pbio.3001007

|

Quick Start
-----------
1\. Install external dependencies

OrthoHMM has two external dependencies — [HMMER](http://hmmer.org/download.html) and [mcl](https://github.com/micans/mcl?tab=readme-ov-file#installation-and-mcl-versions) — that can't be installed using pip.
Download and install these programs from their respective websites, which are linked in the previous sentence.

|

2\. Install OrthoHMM

.. code-block:: shell

	# install
	pip install orthohmm
	# run
	orthohmm <path_to_directory_of_FASTA_files>

Below are more detailed instructions, including alternative installation methods.

|

**1) Installation**

**If you are having trouble installing OrthoHMM, please contact the lead developer, Jacob L. Steenwyk, via [email](https://jlsteenwyk.com/contact.html) or [Bluesky](https://bsky.app/profile/jlsteenwyk.bsky.social) to get help.**

1\. Install external dependencies

OrthoHMM has two external dependencies — [HMMER](http://hmmer.org/download.html) and [mcl](https://github.com/micans/mcl?tab=readme-ov-file#installation-and-mcl-versions) — that can't be installed using pip.
Download and install these programs from their respective websites, which are linked in the previous sentence.

2a\. Install OrthoHMM from pip

To install using *pip*, we recommend building a virtual environment to avoid software dependency issues. To do so, execute the following commands:

.. code-block:: shell

	# create virtual environment
	python -m venv venv
	# activate virtual environment
	source venv/bin/activate
	# install orthohmm
	pip install orthohmm

**Note, the virtual environment must be activated to use orthohmm.**

|

**Install from source**

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the following commands:

.. code-block:: shell

	# download
	git clone https://github.com/JLSteenwyk/orthohmm.git
	cd orthohmm/
	# create virtual environment
	python -m venv venv
	# activate virtual environment
	source venv/bin/activate
	# install
	make install

To deactivate your virtual environment, use the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

Note, the virtual environment must be activated to use orthohmm.

|

2b\. Install OrthoHMM from source

Similarly, to install from source, we recommend using a virtual environment. To do so, use the following commands:

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/orthohmm.git
	cd orthohmm/
	make install

If you run into permission errors when executing *make install*, create a 
virtual environemnt for your installation:

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/orthohmm.git
	cd orthohmm/
	python -m venv venv 
	source venv/bin/activate
	make install

Note, the virtual environment must be activated to use orthohmm.

|

.. **Install from anaconda**

.. To install via anaconda, execute the following command:

.. .. code-block:: shell

.. 	conda install bioconda::orthohmm

.. Visit here for more information: https://anaconda.org/bioconda/orthohmm

.. |

**2) Usage**

To use OrthoHMM in its simpliest form, execute the following command:

.. code-block:: shell

	orthohmm <path_to_directory_of_FASTA_files>

|

^^^^

.. toctree::
	:maxdepth: 4

	about/index
	advanced/index
	performance_assessment/index
	change_log/index
	other_software/index
	frequently_asked_questions/index

^^^^

