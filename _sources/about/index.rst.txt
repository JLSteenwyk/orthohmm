About
=====

^^^^^

**Ortho**\logy inference using **H**\idden **M**\arkov **M**\odels (**OrthoHMM**) was developed as
part of `Jacob L. Steenwyk <https://jlsteenwyk.github.io/>`_'s post-doctoral work. 

Inferring orthology (that is, genes that have shared ancestry) is a major challenge in bioinformatics.
This challenge is exacerbated when sequence divergence is high -- for example, among organisms that
began diverging dozens of millions of years ago -- which makes remote homology difficult to detect.

A popular approach to identify orthology utilizes graph theory. Specifically, an all-to-all comparison
of sequence similarity is conducted using BLAST or DIAMOND and the resulting similarity scores are used
to cluster groups of highly similar sequences together. The resulting groups of sequences are putative
orthologs; thus, these genes are termed orthologous groups of genes (or orthogroups).

Here, we introduce OrthoHMM, which implements Hidden Markov Models to infer orthology. Hidden Markov Models
have a distinct advantage over BLAST and DIAMOND because HMMs incoporate position-specific information
and patterns that enable more sensitivity and specificity for remote homolog detection. Benchmarking
reveals that OrthoHMM outperforms other software for orthology inference.

|

The Developers
--------------

^^^^^

OrthoHMM is developed and maintained by `Jacob L. Steenwyk <https://jlsteenwyk.github.io/>`_
and `Thomas J. Buida III <www.tjbiii.com>`_.

|

|JLSteenwyk|

|GoogleScholarSteenwyk| |GitHubSteenwyk| |BlueSkySteenwyk| |TwitterSteenwyk| 

`Jacob L. Steenwyk <https://jlsteenwyk.github.io/>`_ is a Howard Hughes Medical Institute
awardee of the Life Science Research Foundation at the University of California, Berkeley.
Find out more information at his `personal website <http://jlsteenwyk.github.io/>`_.

.. |JLSteenwyk| image:: ../_static/img/Steenwyk.jpg 
   :width: 35%

.. |GoogleScholarSteenwyk| image:: ../_static/img/GoogleScholar.png
   :target: https://scholar.google.com/citations?user=VXV2j6gAAAAJ&hl=en
   :width: 4.5%

.. |BlueSkySteenwyk| image:: ../_static/img/Bluesky.png
   :target: https://bsky.app/profile/jlsteenwyk.bsky.social
   :width: 4.5%

.. |TwitterSteenwyk| image:: ../_static/img/Twitter.png
   :target: https://twitter.com/jlsteenwyk
   :width: 4.5%

.. |GitHubSteenwyk| image:: ../_static/img/Github.png
   :target: https://github.com/JLSteenwyk
   :width: 4.5%

|

|TJBuida|

|GoogleScholarBuida| |GitHubBuida| |BlueSkyBuida| |TwitterBuida|

`Thomas J. Buida III <http://tjbiii.com/>`_ is a senior software and data engineer at
`Initial State <https://www.initialstate.com/>`_. 
Find out more information at his
`personal website <http://tjbiii.com/>`_.


.. |TJBuida| image:: ../_static/img/Buida.jpeg  
   :width: 35%

.. |GoogleScholarBuida| image:: ../_static/img/GoogleScholar.png
   :target: https://scholar.google.com/citations?user=tXNWGFEAAAAJ&hl=en
   :width: 4.5%

.. |BlueSkyBuida| image:: ../_static/img/Bluesky.png
   :target: https://bsky.app/profile/tjbiii.bsky.social
   :width: 4.5%

.. |TwitterBuida| image:: ../_static/img/Twitter.png
   :target: https://twitter.com/thomasbuida
   :width: 4.5%

.. |GitHubBuida| image:: ../_static/img/Github.png
   :target: https://github.com/TJBIII
   :width: 4.5% 

|

More Team Members
-----------------

^^^^^

|NKing|

|GoogleScholarKing|

`Nicole King <https://kinglab.berkeley.edu/>`_ is a Howard Hughes Medical Institute Investigator
and Professor of Molecular and Cell Biology at the University of California, Berkeley.
Find out more information at her `laboratory’s website <https://kinglab.berkeley.edu/>`_.

.. |NKing| image:: ../_static/img/NKing.jpg
   :width: 35%

.. |GoogleScholarKing| image:: ../_static/img/GoogleScholar.png
   :target: https://scholar.google.com/citations?hl=en&user=PDOSGdIAAAAJ
   :width: 4.5%

|

|ARokas|

|GoogleScholarRokas| |TwitterRokas| 

`Antonis Rokas <https://as.vanderbilt.edu/rokaslab/>`_ is the Cornelius Vanderbilt Chair in 
Biological Sciences and Director of the `Evolutionary Studies Initiative 
<https://www.vanderbilt.edu/evolution/>`_ at `Vanderbilt University <https://www.vanderbilt.edu/>`_.
Research in his laboratory focuses on the study of the DNA record to gain insight into the patterns and 
processes of evolution. Using a combination of computational and experimental approaches, his lab’s current
research aims to understand the molecular foundations of the fungal lifestyle, the reconstruction of the
tree of life, and the evolution of human pregnancy. Find out more information at his 
`laboratory’s website <https://as.vanderbilt.edu/rokaslab/>`_.

.. |ARokas| image:: ../_static/img/Rokas.jpeg
   :width: 35%

.. |GoogleScholarRokas| image:: ../_static/img/GoogleScholar.png
   :target: https://scholar.google.com/citations?user=OvAV_eoAAAAJ&hl=en
   :width: 4.5%

.. |TwitterRokas| image:: ../_static/img/Twitter.png
   :target: https://twitter.com/RokasLab
   :width: 4.5%

|
