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
