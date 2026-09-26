# Installed OrthoBench Readback Preparation

Job 22179 remains a separate, unfinished full inference at the time of this
preparation. The new reader was tested against the historical P1/C1/R1
partition, not unfinished current-run outputs. It parses root-HOG columns
directly rather than extracting gene-like substrings from arbitrary text.
It rejects empty/repeated groups, duplicated membership, foreign genes and
incomplete coverage. Label-invariant comparison reports both unchanged groups
and genes in changed groups, without interpreting a mismatch as a new default.

The [executed baseline readback](installed_ob_baseline_readback_20260926.json)
covers all 251,378 genes and reconstructs the exact retained aggregate score
and all 70 per-RefOG records: F1 74.10607351873405%, precision
81.77045380181866%, recall 67.75533630827974%. This validates parsing/scoring
against the retained evidence, not native inference reproduction or an
independent implementation of the scoring formula. Eleven focused tests pass.

```sh
python -m benchmark_tools.audit_installed_orthobench --repo . \
  --directory benchmarks/work/publication_installed_orthobench_20260926 \
  --baseline-only --output NEW_BASELINE_READBACK.json
python -m pytest -q tests/unit/test_audit_installed_orthobench.py \
  tests/unit/test_run_installed_orthobench.py
```

After successful native/scheduler completion, the `--job 22179` mode checks
the pinned plan and execution report, recorded source/input identities,
coverage, no phylogeny checkpoint reuse, inferred species-tree mode and basic
reconciled-tree file presence. It reports score and partition differences
instead of hiding them. It does not yet parse every gene tree or validate the
pairwise table; `full_phylogeny_validation_complete` remains false. Exact
current input-directory inventory and full output semantics must also be
checked before claiming complete installed-dataset reproduction. No score
admission or historical score replacement follows from this preparation.

The [execution protocol](INSTALLED_ORTHOBENCH_PROTOCOL_20260926.md) remains
unchanged. No job was restarted, no pipeline source/default changed, and no
DGX work or controlled timing measurement occurred.
