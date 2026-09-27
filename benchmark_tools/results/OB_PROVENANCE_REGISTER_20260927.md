# Updated OrthoBench Provenance Register

The new [table](orthobench_provenance_register_20260927/register.md) and
[machine-readable register](orthobench_provenance_register_20260927/register.json)
incorporate the [historical OrthoHMM audit](OB_ORTHOHMM_RETAINED_PROVENANCE_20260927.md).
The earlier September 26 register remains unchanged. All eight scores and
prediction identities, and all six non-OrthoHMM rows, are unchanged by exact
comparison of the two registers.

The added evidence binds each OrthoHMM row to its audited prediction before
adding commands, input evidence and resource metadata. Historical phylogenetic
inference took 3,274.102675 seconds on the shared host. High sensitivity's
319.467197-second cached-hit replay is a separate field, not inserted into
the full-inference column. Initial search time for that row remains unknown.
Process RSS and sampled summed process-tree RSS retain different labels.

The exporter rechecks the hash-pinned audit and its recorded input/output
dependencies before writing. Twenty-four focused register and historical-audit
tests pass, including wrong prediction bindings and prevention of replay-time
substitution or controlled-timing claims. These retained metrics do not imply
complete external-runtime provenance or a controlled speed ranking.

```bash
python -m benchmark_tools.assemble_ob_provenance_register --repo . \
  --output /absolute/new/orthobench-register --include-orthohmm-audit
```

This closes the omission of available historical OrthoHMM metadata from the
eight-method register. It does not supply missing replay FASTA checksums,
unknown comparator resources, comparable hardware measurements or all
cross-dataset provenance. No new inference, score or default selection.
