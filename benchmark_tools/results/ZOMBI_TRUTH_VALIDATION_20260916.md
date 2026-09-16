# Zombi Truth And Extant-Input Validation

`zombi_truth.py` cross-checks root-origin duplication/loss histories without
transfer. It exports only extant proteins, with IDs containing both family
and native lineage identity. Complete simulator FASTAs contain ancestors;
they are not used directly as inference input.

## Cross-Checks

- Parse event TSVs with exact headers, supported event types, finite ordered
  nonnegative times, one root origin, and unique lineage terminal events.
- Parse reconciled XML structurally with ElementTree. Require one tree,
  recognized event types, unique node names, and matching species locations.
- Require event type and child topology to agree between event tables and
  reconciled XML, ignoring biologically irrelevant sibling order.
- Reject disconnected, cyclic, multiply parented, missing, or wrong-arity
  lineages. Speciation descendants cannot overlap in extant species.
- Match surviving terminal events to genes in extant species' native genome
  TSVs and to pruned gene-tree tips. No extant genome gene may be unaccounted.
- Require the extant proteins to be present with unique identifiers and the
  standard amino-acid alphabet in complete sequence output. Extra ancestral
  sequences remain outside inference inputs.

Ortholog pairs are generated only across child branches of speciation
events. Duplication branches do not generate ortholog relations; loss
branches contribute no surviving genes. Tests separately cover duplication
before speciation, duplication after speciation (co-orthologs), and loss.
The `families` field denotes root-origin homolog families, not automatically
root orthogroups: duplications on the species-tree stem can split orthogroups
at the extant species' common ancestor. Do not score those family cliques as
resolved orthology or silently substitute them for a root-HOG reference.

## Validated Smoke Histories

| History | Extant species | Extant genes | Root-origin families | Cross-species ortholog pairs |
| --- | ---: | ---: | ---: | ---: |
| repeat_a, seed 20260916 | 4 | 41 | 10 | 63 |
| independent, seed 20260917 | 4 | 44 | 10 | 66 |

The machine-readable [repeat-a truth](zombi_truth_repeat_a_20260916.json)
and [independent truth](zombi_truth_independent_20260916.json) retain input
hashes, source hash, exact ortholog pairs, family membership, and prepared
FASTA hashes. Raw exports are in `benchmarks/work/zombi_truth_repeat_a_v1`
and `benchmarks/work/zombi_truth_independent_v1`. These results are truth
validation, not OrthoHMM or comparator accuracy scores.

```bash
python benchmark_tools/zombi_truth.py \
  --run benchmarks/work/zombi_smoke_v2/repeat_a \
  --output benchmarks/work/zombi_truth_repeat_a_v1
python benchmark_tools/zombi_truth.py \
  --run benchmarks/work/zombi_smoke_v2/independent \
  --output benchmarks/work/zombi_truth_independent_v1
```

Existing output directories are preserved; use new destinations to rerun.
Fifteen truth-adapter tests include an end-to-end synthetic fixture and
deliberate disagreements in XML, native Newick, genome IDs, and FASTA IDs.
Together with the seed-adapter tests, 20 tests pass. Real extraction passed
on both independent simulator histories.

## Limits And Next Steps

This validates internal consistency and known small truth cases, not the
biological realism of the simulator or every event mode. Transfer and
non-root origination are deliberately rejected. Completely extinct families,
species extinction, and single-tip sequence-output behavior need explicit
fixtures before expanding the scientific panel; missing files must not be
silently treated as absence of truth. XML times are discretized and are not
used as exact event times. Substitution-model rates and branch lengths have
not been independently statistically calibrated here.

Freeze the scientific multi-seed conditions and evaluation units next,
including duplication/loss, divergence, missingness and uneven sampling.
Add deterministic missing-data transformations with their own manifests,
and distinguish those transformations from simulated evolutionary loss.
Run matched methods on extant-only inputs and evaluate event-derived pairs.
Root-HOG/group scoring requires a separately validated root-time definition.
Neither the small smoke panel nor this adapter establishes robustness or
generalization by itself.
