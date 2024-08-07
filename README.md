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
        <a href="https://codecov.io/gh/jlsteenwyk/orthohmm" alt="Coverage">
          <img src="https://codecov.io/gh/jlsteenwyk/orthohmm/branch/master/graph/badge.svg?token=0J49I6441V">
        </a>
        <a href="https://github.com/jlsteenwyk/orthohmm/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/orthohmm">
        </a>
        <a href="https://twitter.com/intent/follow?screen_name=jlsteenwyk" alt="Author Twitter">
            <img src="https://img.shields.io/twitter/follow/jlsteenwyk?style=social&logo=twitter"
                alt="follow on Twitter">
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
        <a href="LINK">
          <img src="https://zenodo.org/badge/DOI/DOI_HERE.svg">  
        </a>   
    </p>
</p>


OrthoHMM infers gene orthology using Hidden Markov Models.<br /><br />
If you found orthohmm useful, please cite *MANUSCRIPT TITLE*. Steenwyk et al. 2020, JOURNAL. doi: [DOI](https://jlsteenwyk.com/publication_pdfs/2020_Steenwyk_etal_PLOS_Biology.pdf).

<br />

---

<br />

This documentation covers downloading and installing OrthoHMM. Details about each function as well as tutorials for using OrthoHMM are available in the [online documentation](https://jlsteenwyk.com/OrthoHMM/).

<br />

**Quick Start**

```shell
# install
pip install orthohmm
# run
orthohmm <path_to_directory_of_FASTA_files>
```

<br />

**Installation**

**If you are having trouble installing OrthoHMM, please contact the lead developer, Jacob L. Steenwyk, via [email](https://jlsteenwyk.com/contact.html) or [twitter](https://twitter.com/jlsteenwyk) to get help.**

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
