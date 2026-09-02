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
search engine that replaces the `phmmer` subprocess. The production CLI was
rebenchmarked with summed process-tree memory on a single 32-core node:

| proteomes | proteins | wall time | peak process-tree RSS | orthogroups |
|----------:|---------:|----------:|----------------------:|------------:|
|         5 |   15,932 |      7.0s |             1.36 GiB |      12,995 |

The output includes singleton orthogroups. The previous 20-100 proteome table
came from a separate experimental driver and was removed pending a production
harness rerun. See `PERFORMANCE_OPTIMIZATION.md` for commands, checksums,
stage timings, and rejected experiments. The legacy `phmmer` path is still available via
`--search_mode phmmer` but is no longer the default.

**Clustering step.** OrthoHMM uses **Leiden CPM with resolution=0.1** as
the default clustering algorithm. This setting is the best local choice for
the highest-priority accuracy benchmarks: it outperforms MCL@1.5 on
OrthoBench and is the recommended fixed setting for QfO-style mixed inputs.
MCL@1.5 remains available as a conservative robustness option for very
distant cross-kingdom datasets, where fixed Leiden CPM can fragment heavily.

Power users can opt into MCL via `--clustering mcl --inflation_value 1.5`,
or use adaptive Leiden via `--clustering leiden --cpm_resolution auto`.

For divergent datasets where recall is more important than runtime, use
`--accuracy_profile high_sensitivity`. This profile broadens the built-in
k-mer candidate search and runs reciprocal-best normalized-hit clustering
followed by singleton reassignment and strict center-star MSA-profile
expansion. Profile candidates must score at least as well as the weakest
existing group member and receive a sequence-supported anchor before they can
change a cluster. The profile currently requires the built-in search and
Leiden clustering. The default `standard` profile is unchanged.

**Phylogeny-aware mode (experimental).** OrthoHMM can reconcile ambiguous
multi-copy families against a validated rooted species tree. This mode is
opt-in; `--phylogeny off` remains the default and does not load tree libraries
or require external phylogenetic programs.

Install the optional parser dependency and provide MAFFT, FastTree, and a
rooted Newick tree whose leaf names match input proteome filenames without
their FASTA extensions:

```shell
pip install 'orthohmm[phylogeny]'
orthohmm proteomes/ --accuracy_profile high_sensitivity \
  --phylogeny reconcile --species_tree_mode supplied \
  --species_tree species_tree.nwk --aligner mafft --tree_builder FastTree
```

To infer the species tree from OrthoHMM's own occupancy-ranked, single-copy
families instead, omit `--species_tree` and use
`--species_tree_mode infer`. Inference uses a high-occupancy marker core and,
when necessary, connected lower-occupancy single-copy markers to place taxa
that are absent from the core. Supplied-tree and inferred-tree results should
be reported separately because they use different external information budgets.

Only multi-copy, multi-species candidate families receive alignments and gene
trees. Results and restart checkpoints are written under
`orthohmm_phylogeny/`, including root-level HOGs, hierarchical node records,
native pairwise orthologs, rooted and annotated gene trees, a reconciliation
summary, and a provenance manifest. The conventional
`orthohmm_orthogroups.txt` remains available and contains the reconciled
root-level groups when this mode is enabled.

Experimental reconciliation policies are independently configurable. The
legacy defaults are `--phylogeny_root_rule supported_children` and
`--phylogeny_pair_rule lca`. The alternative
`--phylogeny_pair_rule positive_paralogy` preserves medium-confidence pairs at
species-tree mapping conflicts and removes pairs only when the two gene-tree
children share an observed species. Pair confidence is written to
`orthohmm_pairwise_orthologs_confidence.tsv`; the original pair file and its
schema are unchanged.
The `--phylogeny_root_rule confidence` policy additionally permits a
root-level split supported by repeated species overlap or a gene-tree branch
support of at least 0.9; missing and weak support remain conservative.
Use `--phylogeny_candidates satellite_v1` with
`--accuracy_profile high_sensitivity` to opt into bounded HMM-backed candidate
expansion. The original groups are retained as seed provenance in
`orthohmm_working_res/phylogeny_candidate_seeds.tsv`.
The `satellite_v2` profile adds seed-bounded iterative expansion and requires
each attachment to receive a high-confidence cross-seed ortholog pair from the
reconciled gene tree before it can remain in the final root HOG. Candidate
evidence is recorded in
`orthohmm_working_res/phylogeny_candidate_merges.json`. A complete v2 run that
infers its own species tree can be requested with:

```shell
orthohmm proteomes/ --accuracy_profile high_sensitivity \
  --phylogeny reconcile --phylogeny_candidates satellite_v2 \
  --species_tree_mode infer --species_tree_rooting min_variance \
  --phylogeny_root_rule species_overlap \
  --phylogeny_pair_rule positive_paralogy \
  --aligner mafft --tree_builder FastTree
```

For internally inferred species trees,
`--species_tree_rooting min_variance` minimizes root-to-tip distance variance
over internal edges and avoids placing the root on a single long terminal
branch when no explicit outgroup is available. The backward-compatible default
remains `midpoint` while this policy is under evaluation.

---

<br />

This documentation covers downloading and installing OrthoHMM. Details about each function as well as tutorials for using OrthoHMM are available in the [online documentation](https://jlsteenwyk.com/orthohmm/).

<br />

**Quick Start**

1\. Install external dependencies

OrthoHMM's default pipeline does not require external bioinformatics
binaries. [MCL](https://github.com/micans/mcl?tab=readme-ov-file#installation-and-mcl-versions)
is optional and only required if you opt into `--clustering mcl`. Install it
via your package manager (`apt install mcl`, `brew install mcl`,
`conda install -c bioconda mcl`) or from source.

[HMMER](http://hmmer.org/download.html) is optional and only required if you opt into the legacy `--search_mode phmmer` pipeline; the default built-in search engine has no HMMER dependency.

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

OrthoHMM's default pipeline does not require external bioinformatics
binaries. [MCL](https://github.com/micans/mcl?tab=readme-ov-file#installation-and-mcl-versions)
is optional and only required if you opt into `--clustering mcl`. Install it
via your package manager (`apt install mcl`, `brew install mcl`,
`conda install -c bioconda mcl`) or from source.

[HMMER](http://hmmer.org/download.html) is optional and only required if you opt into the legacy `--search_mode phmmer` pipeline; the default built-in search engine has no HMMER dependency.

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
