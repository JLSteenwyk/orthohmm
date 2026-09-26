# Installed Frozen-Source Phylogeny Fixture

The [setup-overlay wheel](PUBLICATION_FROZEN_OVERLAY_INSTALL_20260926.md) now
has a separate successful inferred-phylogeny fixture run outside its source
directory. The scientific package remains frozen at `7f3a9e4`, with the same
new compilation and locked dependencies as the preceding installed tests.
No original benchmark executable, scientific default or score was changed.

The [execution report](publication_frozen_phylogeny_smoke_20260926.json) has
SHA-256 `8d258122281bed99d26e4c66bdc03aa3b8bfa3367db3de1bdb5daf86c815e4e0`.
The [separate installed-DendroPy tree readback](publication_frozen_phylogeny_tree_readback_20260926.json)
has SHA-256 `28bf1915e4aeae756255017f1942e88241e39ef025674bff6c14ab100f9633c1`.

## Fixed Fixture and Results

Seed 20260926 generates four taxa, two single-copy marker families and one
family with two identical copies per taxon: 16 proteins, each 200 residues.
These are synthetic integration inputs, not an evolutionary simulation or
biological benchmark. Duplicate sequences intentionally force a multi-copy
candidate while markers permit species-tree inference. The existing 38-gene
fixture lacks single-copy marker coverage for some taxa and was not used to
claim successful species-tree inference.

| Check | Observation |
| --- | ---: |
| Input genes / exactly-once output coverage | 16 / 16 |
| Candidate families / final root HOGs | 3 / 3 |
| Species-tree marker families | 2 |
| Reconciled / bypassed families | 1 / 2 |
| Cross-species ortholog pairs | 36 |
| Reported duplication / speciation events | 4 / 3 |
| Ordinary / remapped gene-tree checkpoint hits | 0 / 0 |
| Species-tree checkpoint hit | false |

The installed CLI ran the built-in HMM high-sensitivity pipeline, satellite_v2
candidate expansion, internally inferred species tree, min_variance rooting,
species_overlap root rule and positive_paralogy pair rule, with one CPU/thread.
It explicitly used local MAFFT 7.525 and FastTree 2.2.0. Both tool entrypoints
are hashed, and the native provenance manifest records resolved paths/versions.
No supplied species tree or pre-existing checkpoint was provided.

The verifier requires nonzero marker and reconciled-family counts, fresh
inference, valid cross-species pairs, complete duplicate-free gene coverage
and nonempty reconciled-tree artifacts. Ten unit tests cover the deterministic
fixture, complete artifacts and rejection of bypass-only/cached/no-marker
results, incomplete or duplicate membership, invalid/duplicate pairs and
missing trees. A separate `python -I` readback using installed DendroPy 5.1.0
parses the rooted species tree and reconciled gene tree: exactly four species
tips and eight expected duplicate-family tips, respectively, without duplicate
leaves. Rooted status is present in both files. This checks tree structure and
tip inventory, not biological correctness of their topology or event labels.

All 159 recorded input/source/installed/output/log record entries were rechecked.
Inputs remained unchanged. The inference log reports 1.545 seconds, a
descriptive shared-host observation, not controlled comparative performance.

## Reproduction

The script uses the hash-pinned preceding installed-source audit and verifies
the isolated package location before execution. It creates a fresh output
directory, runs one inference attempt, limits it to 180 seconds and kills only
that child's process group on timeout. All failures retain a report/log; no
automatic retry occurs.

```sh
/usr/bin/python3 -S -m benchmark_tools.verify_frozen_phylogeny_install \
  --repo . --artifact benchmarks/work/publication_frozen_overlay_20260926 \
  --mafft /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mafft-7.525-with-extensions/bin/mafft \
  --fasttree /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/FastTree_v220/FastTree \
  --output benchmarks/work/publication_frozen_phylogeny_smoke_20260926
python -m pytest -q tests/unit/test_verify_frozen_phylogeny_install.py
```

Choose a new output directory for a separately justified repetition. The
separate tree readback uses `dendropy.Tree.get(..., schema="newick",
preserve_underscores=True)`, checks `is_rooted` and compares leaf-label sets to
`S1` through `S4` and their `family2`/`family2_duplicate` genes. Its raw report
retains tree hashes and the installed DendroPy module path.

## Limits

This executes the full installed pipeline on a small fixture, not the complete
publication datasets, independent biological validation, historical runtime
equivalence or cross-host/platform portability. MAFFT and FastTree remain local
external installations: hashing entrypoints is not a complete transitive
runtime/archive or redistribution review. Identical-copy synthetic data do not
establish difficult paralog resolution. High-CPM remains unadmitted, controlled
timing remains unmet, and the full publication goal is active. No DGX access,
external release or upload occurred.
