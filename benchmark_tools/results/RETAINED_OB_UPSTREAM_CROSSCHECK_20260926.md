# Upstream OrthoBench Scoring Cross-Check

Executed the retained upstream evaluator's input checks and scoring functions
for all five comparator readbacks. The
[report](retained_ob_upstream_crosscheck_20260926.json) retains 108 input/source
records, upstream stdout, exact numerical returns and differences from the
project's official-formula reimplementation. All 15 numerical differences
(five methods times precision, recall and F1) are exactly zero.

The evaluator SHA256 is
`81eb1e660c17819549b07eea8a54b4fb42a89180cafeb4569d92195d282f5e6f`.
Its reference and low-certainty directory inventory must match the earlier
hash-pinned readback, and all 12 input FASTAs are recorded. The expected universe
contains 251,378 genes. Its smart parser is used for the sequence-only
OrthoFinder checkpoint, SonicParanoid, ProteinOrtho and annotated OrthoMCL
output. FastOMA uses the shared two-column adapter before calling upstream
functions; parser independence is not claimed for that format.

Missing-gene messages remain in stdout. FastOMA supplies 130,010 recognized
genes (51.7%); OrthoMCL supplies 216,949 (86.3%). The latter is one fewer than
the 216,950 raw annotated members reported historically. A subsequent
[run inventory comparison](ob_orthomcl_run_gene_inventory_20260926.json) locates
that count difference between the April and July outputs: April includes
`ENSP00000476467`, which July lacks, with no July-only genes. The July annotated
parser and upstream regex both recover 216,949 genes, so a parsing discrepancy
does not explain this count difference. The cross-check here uses July, while
the historical count matches April. Partition equality or score-effect bounds
do not follow from this gene inventory check. SonicParanoid and ProteinOrtho inputs are
singleton-padded scoring files, so their 100% recognized coverage is not native
grouping coverage.

A subsequent [full partition/count comparison](OB_ORTHOMCL_PARTITION_COMPARISON_20260926.md)
finds 23,802 identical groups and one April group split into two July groups
with one missing member. No changed-group genes are in the reference; all
70 family count records and aggregate scores match exactly. Whole-partition
equivalence is explicitly false.

Twenty focused wrapper/scorer/adapter tests pass. A Python SyntaxWarning from
the unchanged upstream regex literal was emitted on first import; the pinned
file was not edited. An inefficient read-only gene-difference probe was stopped
and replaced with a linear set comparison; no inference or scoring run was
restarted. These were function-level upstream executions, not its
CLI wrapper or a new inference run. Historical consumption, native conversion,
timing and complete publication readiness remain unproven.

```bash
/usr/bin/python3 -S -m benchmark_tools.crosscheck_ob_upstream \
  --root . \
  --evaluator /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/Open_Orthobench/BENCHMARKS/benchmark.py \
  --output /tmp/retained-ob-upstream-crosscheck-new.json
```
