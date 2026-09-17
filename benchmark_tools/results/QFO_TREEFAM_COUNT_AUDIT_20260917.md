# QfO TreeFam-A Count and Resampling Audit

The four admitted recovered-stage assessments were checked against the retained
native reference and scorer. This advances scoring verification but does not
complete TreeFam-A uncertainty or the broader QfO ablation requirement.

## Native Reference Semantics

The pinned 2020 reference contains one case, `TreeFamA`, with 11,140 mapped
proteins and 79,320 one-direction relations. Exactly11,130 proteins occur in
those relations; ten mapped proteins have no reference relations. These ten
are not missing predictions and must not be assigned artificial false-negative
counts. The native raw output contains all79,320 relations for each stage.
Reference truth labels and represented member identities agree across stages.

The native scorer adds one uniform-prior count to each confusion category for
the entire pooled case. Its effective counts are `raw_count / 2 + 1`, not
`raw_count + 1`, and not one prior per original TreeFam family. Container and
local scorer hashes match the reviewed source. The independent Darwin read
verifies orientation, endpoint membership counts and absence of reference
endpoints outside the mapped-protein universe.

## Verified Results

Counts below are raw relations before the native prior. P and R are native
precision and recall; F1 is their harmonic mean, not mean per-family F1.

| Recovered stage | TP | FP | FN | TN | P | R | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| multipass | 28,841 | 12,552 | 14,701 | 23,226 | 0.696741 | 0.662357 | 0.679114 |
| multipass_refined | 19,503 | 4,584 | 24,039 | 31,194 | 0.809638 | 0.447917 | 0.576755 |
| strict_profiles | 28,859 | 12,614 | 14,683 | 23,164 | 0.695831 | 0.662770 | 0.678899 |
| strict_profiles_refined | 19,508 | 4,643 | 24,034 | 31,135 | 0.807700 | 0.448032 | 0.576358 |

All reconstructed metrics agree with the admitted native outputs to absolute
tolerance5e-8. The stage name `refined` denotes sequence-based post-clustering
refinement, not phylogenetic reconciliation. Profile contrasts include downstream
responses. These are development-exposed, cluster-derived predictions; the
table does not establish significance, independent generalization or superiority.

## Why No Family Bootstrap Yet

`generateData/AddReconciledTree.drw` iterates through the original `.nhx` trees,
unions their mapped members and merges their relation tables into one
`RecTreeCase`. The raw assessment therefore labels every row `TreeFamA` and
does not retain an original family identifier. The serialized tree argument is
the last loop tree, not a complete substitute for the original tree collection.
No `.nhx` files or `treefam2reference` mapping were found by filename inventory
under the local `qfo_benchmark` directory. This is a local availability finding,
not proof that the upstream data cannot be recovered.

Resampling the one pooled case is degenerate. Resampling relations independently
ignores their shared proteins and evolutionary history. Reference-graph connected
components have not been validated as original independent families. None of
these shortcuts is reported as a family-level confidence interval.

Next recover the original family-to-reference mapping with source provenance,
verify how overlaps and overwritten relations were handled, and freeze a
resampling protocol before inspecting family-level contrasts. Any resampled
statistic must reconstruct the pooled native endpoint with its one prior, not
substitute a macro-family statistic. SwissTrees intervals do not transfer to
this challenge; the six-metric secondary mean also remains without a justified
joint uncertainty model.

## Evidence and Execution

[Machine-readable audit](qfo_treefam_counts_20260917.json) records all raw counts,
source hashes, container identity and Darwin inventory output. Reference SHA256:
`6f419f96886e5cf14ac889a1437cdb2655140cb8fd281f6bb5f8e0baa69d23b3`.
The initial equality check between mapped and raw-represented proteins failed
and prompted the isolated-member check. A first endpoint-set implementation
timed out after60s; using a table to accumulate unique endpoints completed.
Neither failure changed inference, reference truth or native scores.

Thirteen focused count-reader/metric tests pass. Current raw-file hashes and
generator source are retrospective provenance, not execution-time attestations.
Run from the repository with retained native artifacts and a fresh output:

```sh
python benchmark_tools/audit_qfo_treefam_counts.py --repo . --output /tmp/qfo_treefam_counts.json
```
