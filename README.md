<p align="center">
  <img src="https://raw.githubusercontent.com/JLSteenwyk/orthohmm/main/docs/_static/img/logo.png" alt="OrthoHMM" width="400">
</p>

# OrthoHMM

HMM-centered orthology inference from protein FASTA files.

[Documentation](https://jlsteenwyk.com/orthohmm/) |
[Report an issue](https://github.com/JLSteenwyk/orthohmm/issues) |
[CI](https://github.com/JLSteenwyk/orthohmm/actions) |
[Coverage](https://codecov.io/gh/JLSteenwyk/orthohmm)

OrthoHMM combines a built-in k-mer prefilter and profile-HMM scoring with
Leiden clustering. High-sensitivity inference and phylogenetic reconciliation
are optional. The default built-in search does not require HMMER or MCL.

## Install From Source

Python 3.10 or newer is required. A C compiler with OpenMP support builds the
native CPU kernels during installation. Without it, standard search can use
the slower Numba fallback; high-sensitivity profile expansion requires the
native alignment kernel. Linux users can install GCC and its development
runtime through their system package manager. On macOS, GNU GCC/OpenMP must
be configured before installation for native kernels.

```sh
git clone --depth 1 https://github.com/JLSteenwyk/orthohmm.git
cd orthohmm
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
orthohmm --help
```

This installs the checked-out revision, which may differ from the published
PyPI release. For a generic-CPU source build rather than a host-optimized
build, use `ORTHOHMM_CPU_TARGET=baseline python -m pip install .`.
Built wheels remain platform-specific. CUDA compilation is optional when
`nvcc` is available; GPU hardware is not required for CPU inference.

## Run

Place one protein FASTA file per species in a directory. Accepted extensions
are `.fa`, `.faa`, `.fas`, `.fasta`, `.pep`, and `.prot`. Sequence identifiers
must be unique within each file. Create the output directory before running.

```sh
mkdir results
orthohmm proteomes/ -o results/ -c 8
```

Outputs include `orthohmm_orthogroups.txt`, per-group FASTA files, gene counts,
and working files. Input FASTA files are not modified. To increase search
sensitivity at additional computational cost:

```sh
mkdir results-sensitive
orthohmm proteomes/ -o results-sensitive/ -c 8 --accuracy_profile high_sensitivity
```

MCL clustering and legacy sequence search remain available when their external
programs are installed:

```sh
mkdir results-mcl results-phmmer
orthohmm proteomes/ -o results-mcl/ --clustering mcl --inflation_value 1.5
orthohmm proteomes/ -o results-phmmer/ --search_mode phmmer
```

## Optional Phylogeny

Install the tree-parser dependency and provide MAFFT and FastTree on PATH.
Species-tree leaf names must match proteome filenames without FASTA extensions.

```sh
python -m pip install '.[phylogeny]'
mkdir results-phylogeny
orthohmm proteomes/ -o results-phylogeny/ -c 8 \
  --accuracy_profile high_sensitivity --phylogeny reconcile \
  --species_tree_mode supplied --species_tree species_tree.nwk \
  --aligner mafft --tree_builder FastTree
```

To infer the species tree instead, use `--species_tree_mode infer` and omit
`--species_tree`. Supplied and inferred trees use different information and
should be reported separately. Phylogenetic inference is not enabled by default.

Optional HMMER profile-control Python helpers require
`python -m pip install '.[hmmer]'` and the corresponding HMMER executables.
See the [advanced documentation](https://jlsteenwyk.com/orthohmm/advanced/)
and `orthohmm --help` for supported options and experimental-mode restrictions.

## Development

```sh
python -m pip install . -r tests/requirements.txt
make kernels
make test
make test.coverage
```

The [test guide](tests/README.md) explains optional external tools, native
kernel preparation, and coverage. CI checks Linux/macOS, Python 3.10-3.13,
source installation, native CPU kernels, and built-wheel inference.

## Research Materials

The application branch does not include manuscripts, benchmark datasets,
analysis outputs, or historical reproduction workflows. Those materials and
their tests remain on
[research-archive-20261029](https://github.com/JLSteenwyk/orthohmm/tree/research-archive-20261029),
preserved at commit `ebfa635ca786796287c9ff2a3caba1c2618f016e`.
The repository's published history has not been rewritten. A normal clone also
checks out only the application, but downloads historical objects; the shallow
clone above avoids downloading that research history.

## Citation And License

Please cite the original method: Steenwyk et al. (2024), *OrthoHMM: Improved
Inference of Ortholog Groups using Hidden Markov Models*, bioRxiv,
[doi:10.1101/2024.12.07.627370](https://www.biorxiv.org/content/10.1101/2024.12.07.627370).
Machine-readable metadata is in [CITATION.cff](CITATION.cff). Report the exact
software revision and inference settings; a paper citation does not identify
a development configuration or establish superiority on arbitrary datasets.

OrthoHMM is distributed under the [MIT license](LICENSE.md).
