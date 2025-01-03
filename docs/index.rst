.. image:: _static/img/logo_condensed.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/orthohmm

^^^^^


OrthoHMM using high sensitivity and specificity Hidden Markov Models for orthology inference.

If you found OrthoHMM useful, please cite *OrthoHMM: Improved Inference of Ortholog Groups using Hidden Markov Models*. Steenwyk et al. 2024, bioRxiv. doi: |doiLink|_.

.. _doiLink: https://www.biorxiv.org/content/10.1101/2024.12.07.627370v1
.. |doiLink| replace:: 10.1101/2024.12.07.627370

|

Quick Start
-----------
1\. Install external dependencies

OrthoHMM has two external dependencies — |hmmerLink|_ and |mclLink|_ — that can't be installed using pip.
Download and install these programs from their respective websites, which are linked in the previous sentence.

.. _hmmerLink: http://hmmer.org/download.html
.. |hmmerLink| replace:: HMMER
.. _mclLink: https://github.com/micans/mcl?tab=readme-ov-file#installation-and-mcl-versions
.. |mclLink| replace:: mcl

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

**If you are having trouble installing OrthoHMM, please contact the lead developer, Jacob L. Steenwyk, via |contactSteenwyk|_ or |blueskySteenwyk|_ to get help.**

.. _contactSteenwyk: https://jlsteenwyk.com/contact.html
.. |contactSteenwyk| replace:: email
.. _blueskySteenwyk: https://bsky.app/profile/jlsteenwyk.bsky.social
.. |blueskySteenwyk| replace:: Bluesky

1\. Install external dependencies

OrthoHMM has two external dependencies — |hmmerLink|_ and |mclLink|_ — that can't be installed using pip.
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

