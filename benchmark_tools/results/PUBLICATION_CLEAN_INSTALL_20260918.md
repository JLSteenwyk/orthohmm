# Clean Dependency Installation Check

Built a wheel from a fresh `git archive 4fe2c28` source tree, then installed
`orthohmm[phylogeny]` into a new Python3.10.13 venv without system site
packages. No inference environment was upgraded. All installed distributions
and the imported OrthoHMM module resolve inside this venv. `pip check`
reported no broken requirements.

Evidence directory: `benchmarks/work/publication_clean_install_4fe2c28/`.
[Machine-readable record](publication_clean_install_4fe2c28_20260918.json)
contains exact resolved versions, wheel/input/output/log hashes and limits.
The wheel SHA-256 is
`b06264e52edb784a6524a1a3a97c2099dbbd2e3a52cb5cbe57ac32d56b059e5b`.

## Commands And Results

From the evidence directory after creating the source archive:

```sh
python3 -m pip wheel --no-deps --no-build-isolation ./source --wheel-dir ./wheels --log ./wheel-build.log
python3 -m venv venv
./venv/bin/python -m pip install 'wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl[phylogeny]' --log ./clean-install.log
./venv/bin/python -m pip check
./venv/bin/orthohmm inputs -o outputs_standard -c 1 --search_mode builtin --clustering leiden
./venv/bin/orthohmm inputs -o outputs_high -c 1 --search_mode builtin --clustering leiden --accuracy_profile high_sensitivity
```

Inputs are copies of the11top-level sample FASTAs from the clean source
archive; output directories were created empty. Both installed CLI runs
exited0 and produced4groups, each covering all38input genes exactly once.
Their partition bytes match on this fixture. An isolated-interpreter check
verified venv-only distributions and used the installed FASTA parser to
reconstruct the input universe. The optional DendroPy dependency also parsed
a two-leaf Newick fixture. This is not a complete phylogenetic inference run.

Resolved runtime versions include NumPy2.2.6, Numba0.67.0, llvmlite0.49.0,
python-igraph/igraph1.0.0, leidenalg0.12.0, texttable1.7.0 and DendroPy5.1.0.
The record is an observed version inventory, not a hash-complete portable
lock; rerunning the unconstrained installation later may select new versions.

## Remaining Limits

The package build used the host's existing build backend and compiler, so
only dependency installation was isolated. Native library/ISA portability,
other Python/OS combinations, a fully pinned fresh build, and external
MAFFT/FastTree execution remain separate release requirements. These small
fixture outputs do not validate biological accuracy or establish equivalence
between dependency versions on publication datasets. No wheel was uploaded
or presented as a completed release; timing is not scientific evidence.
