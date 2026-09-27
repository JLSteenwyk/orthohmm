# Recovery Installation Phylogeny Fixture

The fresh hash-locked recovery venv completed inferred species/gene phylogeny
on the retained 16-gene, four-species fixture. The
[result receipt](publication_recovery_phylogeny_20260926.json) binds the new
plan, exact interpreter/command, native completion and independent readbacks.
All alignments and trees were newly generated. Species-tree checkpoint reuse
is false and both family checkpoint counters are zero.

The result contains three root groups and 36 ortholog pairs: one reconciled
family, two bypassed families, four duplication events, three speciation
events and zero uncertain events. Two families contributed to species-tree
inference. Structural, sequence-content, reconciliation/pair and hierarchy
readers all passed. The complete installation audits before and after native
execution have the same SHA256 as the earlier version-2 installation audit:
`f1636e6e3222c135d2f174bb3f2c247f74e646d8653d087f068960f4926f2559`.
The focused fixture, installation and native-admission suite passes 39 tests.

## Executed Workflow

```bash
python -m benchmark_tools.verify_recovery_phylogeny --repo . \
  --directory /absolute/fresh/recovery-phylogeny
```

The executed output is `benchmarks/work/publication_recovery_phylogeny_20260926`.
Existing output directories are rejected. The runner checks the pinned
historical fixture plan, its recorded inputs/external executables and the
fresh package installation, then substitutes only the new interpreter/runtime
and output directory. Frozen inference settings and scientific code remain
unchanged. Imports must reside inside the fresh installation. The native
launcher permits one attempt and no retry/resume; all native logs are retained.

This is a same-host fixture validation, not a full-data benchmark, biological
accuracy evaluation, complete search-to-phylogeny run or portable installer.
It depends on the retained local template and external MAFFT/FastTree tools.
The previously recovered full OrthoBench score is not transferred to this
environment. Integration of canonical ordering into an executable full-pipeline
reproduction path and validation of that path remain open.
