# Simulation Transform And Scorer Checks

The frozen selection rules in
`PUBLICATION_SIMULATION_PROTOCOL_20260916.md` are implemented in
`simulation_conditions.py`. Selection accepts only gene IDs, their species,
the species tree, the master seed and condition. It cannot inspect family
labels, true pairs, sequence similarity or predictions.

`derive_simulation_conditions.py` first validates baseline native truth,
selects removals, projects true pairs to retained genes, then writes new
extant-only FASTAs and per-condition truth files. It records all removed
genes/species, parent truth/input hashes, transformation source hashes and
export hashes. Existing outputs are rejected. An inapplicable clade rule
is reported explicitly rather than replaced with another selection rule.
No supplied tree is exported; primary methods infer trees from each input.

## End-To-End Smoke Check

Applying the transformations to the four-species seed-20260916 smoke input
produced these audited exports:

| Condition | Retained genes | Eligible true ortholog pairs |
| --- | ---: | ---: |
| baseline smoke input | 41 | 63 |
| missing20 | 33 | 38 |
| uneven_taxa | 31 | 32 |
| taxon_count_control | 30 | 30 |

Missing20 removes exactly floor(41/5) = 8 genes. Uneven sampling and its
count control remove the same number of taxa but need not retain the same
number of genes or eligible pairs. This four-species exercise is a workflow
test, not a member of the frozen eight-species scientific panel. It has no
method accuracy results.

[Transformation manifest](zombi_transforms_smoke_20260916.json) retains
selections, baseline evidence and output truth hashes. Raw exports are
`benchmarks/work/zombi_transforms_smoke_v1`. Reproduction:

```bash
python benchmark_tools/derive_simulation_conditions.py \
  --baseline-run benchmarks/work/zombi_smoke_v2/repeat_a \
  --output benchmarks/work/zombi_transforms_smoke_v1 --seed 20260916
```

Use a new output path when rerunning; previous evidence is preserved.

## Scoring Invariants

Pair scoring includes every predicted cross-species pair in the retained
input universe, including predictions between unrelated origin families.
An adversarial fixture merging two unrelated two-species families yields
two true positives and two false positives, not perfect precision.
Unknown endpoints, same-species/self pairs, duplicate group membership,
and duplicate truth pairs are rejected. Prediction orientation and duplicate
rows are canonicalized, with the number removed recorded explicitly.

Coverage here is the fraction of input genes in at least one predicted pair;
it is not final-group membership coverage or proof that a singleton was
unprocessed. Undefined ratios are represented as zero with explicit metric
names in `undefined_ratios`, not silently treated as measured perfect
performance. The caller must establish successful method completion before
scoring: this function cannot distinguish an empty successful prediction
from an absent failed run. Failed runs must not be replaced by empty pairs.

Thirty focused tests pass across transformation, export, scoring and native
truth checks. Tests cover order-independent hashing, clade tie-breaking,
count matching, inapplicable cases, truth projection, cross-family false
positives, pair orientation, export hashes and non-overwriting behavior.

The executable scientific launch manifest, full panel generation, seed-level
aggregation, native method conversion and real method runs remain open.
The smoke checks do not establish scientific robustness or satisfy the
publication goal by themselves.
