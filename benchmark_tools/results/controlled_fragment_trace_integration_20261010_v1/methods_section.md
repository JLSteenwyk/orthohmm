### Retained Fragment Representative Protocol

Representative selection is retrospective: aggregate results had already been
inspected. Four methods, three truncated-endpoint strata and seven categories
define 84 bins. Categories are fragment FN/FP, new FN/FP relative to baseline,
recovered FN, removed FP and retained TP controls. Within each nonempty bin,
the minimum SHA256 of UTF-8 `method:category:seed:gene_a:gene_b`, with lexical
tie-breaking, selects one eligible pair-record across the ten retained seeds.
Thirty-one bins are empty, without replacement. Fifty-three selected records
are unique in this realization; the protocol permits category overlap and
does not make these observations independent or representative of prevalence.

Read-only tracing retains both arms for every selected case. Significant-hit
NumPy checkpoints provide exact directed rows/scores; the weighted graph and
groups provide direct edges, connectivity and seed membership. PHY additionally
uses saved candidate families, observed reconciliation-node ancestors, logged
merge constraints and final RootHOGs. Positive-paralogy excludes observed
duplications but permits mapping-only uncertain calls; membership filtering
is active only when saved constraints are nonempty. Final RootHOG separation
alone is not proof of native pair exclusion. Unambiguous families without tree
inference are separately marked as bypasses. OF is traced only through its
MCL checkpoint and native pair output, not a matched search or event adapter.

Independent native-file readers reproduce all 106 observations in 44 contexts,
including graph connectivity by BFS and bottom-up root/constraint reconstruction.
The initial readback refused the TSV's null serialization after native checks;
a separate tested correction accepts the original `NA` cells, not empty or
false fields, while preserving the failed reader and original report/table.
No inference, scoring, bootstrap or tuning is repeated.

