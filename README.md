<p align="center">
  <a href="https://github.com/jlsteenwyk/orthohmm">
    <img src="https://raw.githubusercontent.com/JLSteenwyk/orthohmm/master/docs/_static/img/logo.png" alt="Logo" width="400">
  </a>
  <p align="center">
    <a href="https://jlsteenwyk.com/orthohmm/">Docs</a>
    ·
    <a href="https://github.com/jlsteenwyk/orthohmm/issues">Report Bug</a>
    ·
    <a href="https://github.com/jlsteenwyk/orthohmm/issues">Request Feature</a>
  </p>
    <p align="center">
        <a href="https://github.com/JLSteenwyk/orthohmm/actions" alt="Build">
            <img src="https://img.shields.io/github/actions/workflow/status/JLSteenwyk/orthohmm/ci.yml?branch=main">
        </a>
        <a href="https://codecov.io/gh/JLSteenwyk/orthohmm" >
          <img src="https://codecov.io/gh/JLSteenwyk/orthohmm/graph/badge.svg?token=YEXCJN8D4E"/>
        </a>
        <a href="https://github.com/jlsteenwyk/orthohmm/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/orthohmm">
        </a>
        <a href="https://bsky.app/profile/jlsteenwyk.bsky.social" target="_blank" rel="noopener noreferrer">
          <img src="https://img.shields.io/badge/Bluesky-0285FF?logo=bluesky&logoColor=fff">
        </a>
        <br />
        <a href="https://pepy.tech/badge/orthohmm">
          <img src="https://static.pepy.tech/personalized-badge/orthohmm?period=total&units=international_system&left_color=grey&right_color=blue&left_text=PyPi%20Downloads">
        </a>
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <a href="https://pypi.org/project/orthohmm/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/orthohmm">
        </a>
        <a href="https://www.biorxiv.org/content/10.1101/2024.12.07.627370">
          <img src="https://zenodo.org/badge/DOI/10.1101/2024.12.07.627370.svg">  
        </a>   
    </p>
</p>


OrthoHMM infers gene orthology using Hidden Markov Models.<br /><br />
If you found orthohmm useful, please cite *OrthoHMM: Improved Inference of Ortholog Groups using Hidden Markov Models*. Steenwyk et al. 2024, bioRxiv. doi: [10.1101/2024.12.07.627370](https://www.biorxiv.org/content/10.1101/2024.12.07.627370v1).

---

**Performance**

As of v0.2.0, OrthoHMM ships a built-in profile HMM + k-mer prefilter
search engine that replaces the `phmmer` subprocess. It scales to 100
bacterial proteomes (~352K total proteins) on a single 32-core node:

| proteomes | wall time | peak RAM | orthogroups |
|----------:|----------:|---------:|------------:|
|         5 |       13s |  0.29 GB |       1,196 |
|        20 |     4 min |  0.44 GB |       8,680 |
|        60 |    28 min |  1.65 GB |      19,029 |
|       100 |    77 min |  4.67 GB |      27,328 |

Numbers from the bacterial scaling benchmark (RefSeq, 32 threads, 1 TB
node). The legacy `phmmer` path is still available via
`--search_mode phmmer` but is no longer the default.

v0.3.0 also swaps **MCL → Leiden CPM** for the clustering step. On
identical RBNH edges from OrthoBench 2020, Leiden CPM (resolution=0.1)
hits **F=65.7%** vs MCL's best **F=62.4%** at inflation=1.5. MCL remains
available via `--clustering mcl` for users who need bit-for-bit parity
with prior releases.

---

<br />

This documentation covers downloading and installing OrthoHMM. Details about each function as well as tutorials for using OrthoHMM are available in the [online documentation](https://jlsteenwyk.com/orthohmm/).

<br />

**Quick Start**

1\. Install external dependencies

As of v0.3.0, OrthoHMM has **no required external binaries** — clustering uses Leiden via `igraph`/`leidenalg` (installed by pip), and the search engine is built in.

[HMMER](http://hmmer.org/download.html) is optional and only required if you opt into the legacy `--search_mode phmmer` pipeline; [mcl](https://github.com/micans/mcl?tab=readme-ov-file#installation-and-mcl-versions) is optional and only required if you opt into `--clustering mcl`.

<br>

2\. Install OrthoHMM

```shell
# install
pip install orthohmm 
# run
orthohmm <path_to_directory_of_FASTA_files>
```

<br />

**Installation**

**If you are having trouble installing OrthoHMM, please contact the lead developer, Jacob L. Steenwyk, via [email](https://jlsteenwyk.com/contact.html) or [Bluesky](https://bsky.app/profile/jlsteenwyk.bsky.social) to get help.**

1\. Install external dependencies

As of v0.3.0, OrthoHMM has **no required external binaries** — clustering uses Leiden via `igraph`/`leidenalg` (installed by pip), and the search engine is built in.

[HMMER](http://hmmer.org/download.html) is optional and only required if you opt into the legacy `--search_mode phmmer` pipeline; [mcl](https://github.com/micans/mcl?tab=readme-ov-file#installation-and-mcl-versions) is optional and only required if you opt into `--clustering mcl`.

<br>

2a\. Install OrthoHMM from pip

To install using *pip*, we recommend building a virtual environment to avoid software dependency issues. To do so, execute the following commands:
```shell
# create virtual environment
python -m venv venv
# activate virtual environment
source venv/bin/activate
# install orthohmm
pip install orthohmm
```
**Note, the virtual environment must be activated to use *orthohmm*.**

After using OrthoHMM, you may wish to deactivate your virtual environment and can do so using the following command:
```shell
# deactivate virtual environment
deactivate
```

<br />

2b\. Install OrthoHMM from source

Similarly, to install from source, we recommend using a virtual environment. To do so, use the following commands:
```shell
# download
git clone https://github.com/JLSteenwyk/orthohmm.git
cd orthohmm/
# create virtual environment
python -m venv venv
# activate virtual environment
source venv/bin/activate
# install
make install
```
To deactivate your virtual environment, use the following command:
```shell
# deactivate virtual environment
deactivate
```
**Note, the virtual environment must be activated to use *orthohmm*.**

<!-- <br />

To install via anaconda, execute the following command:

``` shell
conda install bioconda::orthohmm
```
Visit here for more information: https://anaconda.org/bioconda/orthohmm -->
