# User-Facing Documentation Reconciliation

30 September 2026. Documentation only; no inference, scoring, dependency
installation, scientific default change or production timing run.

The user cannot currently identify a quiet Threadripper window. Timing is
deferred; other publication work continues. No unrelated process or service
was stopped, and no DGX operation was performed.

## Corrections

- README and Sphinx index distinguish the published release from development
  features and the frozen scientific method. The experimental phylogeny
  example now installs `.[phylogeny]` from an isolated source checkout.
- Sphinx index replaces the removed experimental 20-100-proteome scaling
  table with the historical production measurement already in README and
  `PERFORMANCE_OPTIMIZATION.md`: five proteomes, 15,932 proteins, about 7.0 s,
  1.36 GiB sampled summed process-tree RSS, 12,995 groups including singletons.
  This is not a controlled comparison with competitors or a new observation.
- Removed unqualified optimal-clustering and cross-kingdom accuracy advice.
  Default resolution 0.1 remains unchanged and development-exposed.
- The performance page no longer says dedicated-machine timing is underway.
  README and documentation route readers to corrected QfO version 7 instead
  of the historical 16 September table. Historical reports are unchanged.
- Current benchmark limitations, output semantics, reproduction instructions
  and incomplete release/deposition status are visible from the project entrypoints.

## Published Release Check

Read public primary metadata from
`https://pypi.org/pypi/orthohmm/json` and
`https://pypi.org/pypi/orthohmm/0.2.0/json` on this date. Latest version 0.2.0,
no declared extras (`provides_extra` is null). Runtime requirements are
`numpy>=1.23`, `numba>=0.58`, `python-igraph`, and `leidenalg`.
Published files were uploaded 26 April 2026.

Downloaded the official 438,996-byte `orthohmm-0.2.0.tar.gz` into memory and
checked its SHA256 against the PyPI metadata:
`448a74ec524bb6d1f91cf5ecea13da9d84a0218d1a5b1dafea8f11f908dadaa1`.
Read `orthohmm/parser.py` and `orthohmm/args_processing.py` directly from the
archive without extracting, importing or installing it. Its parser lacks
`--accuracy_profile` and `--phylogeny`; its argument processing requires MCL.
This supports the release distinction, not a complete published-source audit.
The current checkout is labeled 0.5.0 and declares a DendroPy phylogeny extra;
no new distribution was uploaded or release created.

## Validation

The existing documentation built successfully before edits with Sphinx 8.1.3.
The edited documentation built successfully with warnings treated as errors:

```sh
/home/bizon/anaconda3/bin/python -m sphinx -W --keep-going -b html \
  docs benchmarks/work/publication_docs_reconciliation_20260930/after_html
```

All seven RST pages built without warnings. `git diff --check` passed for the
edited pages. Current parser source confirms unchanged defaults: Leiden,
CPM resolution 0.1, standard accuracy profile and phylogeny off. No native
scientific run or installation validation was repeated for these text edits.
Generated HTML remains a local build, not an uploaded website.

Repository documentation changes do not update the existing PyPI description
or prove that the public documentation site has been redeployed. Public release,
archival deposition, complete uncertainty/rights and controlled resource
evidence remain open. The full publication goal is not complete.
