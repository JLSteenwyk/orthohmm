# Current Cross-Dataset Score Table

[Readable scores](current_benchmark_scores_20260926/scores.md),
[exact TSV](current_benchmark_scores_20260926/scores.tsv), and
[source/semantics manifest](current_benchmark_scores_20260926/manifest.json)
consolidate eight methods across OrthoBench, corrected QfO and Three Kingdoms.
The exporter uses three hash-pinned retained reports; it does not rescore data.
The original-release QfO fields in the historical OrthoBench source are never
selected. The obsolete SonicParanoid Three Kingdoms diagnostic is excluded in
favor of its contemporary matched-input result.

All values are shown in 0-to-1 units. OrthoBench weighted F1 is converted from
percent. QfO's GO/EC/FAS similarities remain distinct from its three F1 endpoints;
its six-endpoint mean is a project-defined secondary summary. Three Kingdoms
is a supplementary BUSCO-reference pair test, not genome-wide orthology truth.
Do not average across columns or interpret the table as a universal ranking.
Per-method prediction semantics and Three Kingdoms input-consumption limitations
are retained in the manifest.

Nine tests pass, covering exact source values/units, corrected-only selection,
contemporary SonicParanoid selection, duplicate/missing/admission/metric errors,
source-digest rejection and output preservation. Initial execution exposed that
five historical OrthoBench rows lack direct prediction-file hashes in their
summary objects. The exporter now retains each complete source row and documents
that provenance gap; it does not invent hashes or call this transitive auditing.

A subsequent [five-comparator prediction readback](RETAINED_OB_COMPARATOR_READBACK_20260926.md)
now supplies current file hashes and reproduces those retained statistics.
It does not establish historical consumption or complete native provenance.

```bash
/usr/bin/python3 -S -m benchmark_tools.export_current_benchmark_scores \
  --results benchmark_tools/results --output /tmp/current-benchmark-scores-new
python -m pytest -q tests/unit/test_export_current_benchmark_scores.py
```

No scientific defaults, predictions, uncertainty, timing or native admissions
changed. Full command/version/resource and transitive-provenance consolidation,
controlled timing and the remaining publication requirements are unfinished.
