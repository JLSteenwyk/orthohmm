# Canonical Phylogeny Admission Gate

Prepared while job 22324 was running, without reading its final predictions or
changing the frozen inference code, plan, runtime or input files.

`benchmark_tools.admit_canonical_ob_phylogeny` requires successful terminal
Slurm accounting and a matching submission before inspecting full-run results.
It verifies the frozen plan, all checked input/source records, exact FASTA
directory inventory, start/runtime record, exact native/replay commands,
completion receipt and hashed output bindings. It rejects changed reconciliation
rules, supplied species trees, unconstrained inference and checkpoint reuse.
These are provenance/scope gates, not proof of biological correctness.

The retained 16-gene fixture passes this artifact gate. Its
[receipt](canonical_ob_phylogeny_fixture_admission_20260926.json) deliberately
sets `scientific_scores_admitted` to false. The separate four-reader fixture
validation is recorded in the original
[preparation receipt](canonical_ob_phylogeny_preparation_20260926.json).

112 focused tests pass across admission, strict input/partition readers,
launcher and the four independent phylogeny readers. Mutation tests reject
changed artifacts, nonterminal/failed scheduler states, and semantic changes
even when their output hashes are consistently regenerated. A synthetic test
also explicitly proves that successful artifact binding alone does not validate
group contents or admit scores.

## Complete Post-Run Workflow

`benchmark_tools.readback_canonical_ob_phylogeny` now chains native admission,
all four independent readers, strict full-coverage root-group loading and
three-way frozen scoring. It recomputes both retained baseline scores and
requires exact agreement before comparing canonical results. It preserves
partition changes and every changed reference-family record even when overall
scores agree. It refuses an existing output directory and retains partial
audit files if a later check fails.

```bash
python -m benchmark_tools.readback_canonical_ob_phylogeny \
  --repo . \
  --directory benchmarks/work/canonical_ob_phylogeny_20260926 \
  --job 22324 \
  --output benchmarks/work/canonical_ob_phylogeny_20260926/readback
```

Run only after successful native completion. This command includes admission;
the individual admission command below is an alternative preliminary check,
not a requirement to rerun inference. A comparison preflight recomputed both
retained full-data score objects exactly: historical F1 74.10607351873405% and
fresh-installed F1 73.82156906618316%, with all 70 RefOG records matching.
See the [baseline preflight receipt](canonical_ob_readback_baseline_preflight_20260926.json).
No canonical final result was inspected during this preflight.

The expanded suite passes 115 tests. Its initial run had three failures because
new comparison fixtures used mutable sets whereas the strict partition reader
returns frozensets. The fixtures were corrected to the existing contract; the
scientific scorer and partition-comparison implementation were not changed.

After job 22324 completes, run from the repository root:

```bash
python -m benchmark_tools.admit_canonical_ob_phylogeny \
  --directory benchmarks/work/canonical_ob_phylogeny_20260926 \
  --job 22324 \
  --output benchmark_tools/results/canonical_ob_phylogeny_admission_22324.json
```

Then run the independent structure, sequence, event/pair and hierarchy readers
with the plan's exact input and constraint paths. Only after these checks may
the frozen root partition scorer compare the 70 RefOGs against historical and
fresh-installed runs. Keep all discrepancies and historical results. Do not
interpret this partial-stage diagnostic as a new independent generalization
test, full end-to-end timing run or production default change.
