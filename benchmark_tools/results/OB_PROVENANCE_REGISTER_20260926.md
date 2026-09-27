# OrthoBench Provenance Register

The [readable register](orthobench_provenance_register_20260926/register.md)
and [machine-readable register](orthobench_provenance_register_20260926/register.json)
join all eight historical score rows to six hash-pinned retained reports.
Prediction files are checked against their recorded identities. Scores are
copied exactly from the frozen cross-dataset manifest, not recomputed.

The register consolidates output semantics, available native coverage,
singleton padding, commands/settings, input identity and resource evidence.
It is explicitly partial: historical OrthoHMM full-run provenance/resources
are not filled from the fresh installed run or cached replay measurements.
Unknown values remain null, not zero. The six source reports and the assembler
are recorded with file hashes; this does not attest the full historical
execution chain or historical executable identities.

Proteinortho and the scored July OrthoMCL run each delete 869 asterisks from
177 proteins. These rows cannot be called exact-input comparisons. OrthoMCL's
prediction identity selects July 2025, with eight BLAST threads and a
371,702-second native-log interval; April's matching inputs and 32-thread
record must not be borrowed. The effect of preprocessing on accuracy is unknown.

OrthoFinder checkpoint conversion took 0.42 seconds, but its sequence-only
inference duration remains unknown. FastOMA's 3,023-second workflow includes
failed attempts and used a supplied tree. Tool-reported durations, native-log
intervals and GNU-time measurements have different scopes. None establishes
a controlled comparative runtime ranking; peak process RSS is not aggregate
concurrent-process memory.

Reproduce into a new directory without replacing retained evidence:

```bash
python -m benchmark_tools.assemble_ob_provenance_register \
  --repo . --output /tmp/orthobench-provenance-register
python -m pytest -q tests/unit/test_assemble_ob_provenance_register.py
```

Full cross-dataset provenance consolidation, controlled resource evidence
and the remaining publication requirements are unfinished. No defaults,
predictions or historical benchmark scores changed.
