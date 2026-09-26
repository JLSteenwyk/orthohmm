# Current Cross-Dataset Score Table

[Readable scores](current_benchmark_scores_20260926_v2/scores.md),
[exact TSV](current_benchmark_scores_20260926_v2/scores.tsv), and
[source/semantics manifest](current_benchmark_scores_20260926_v2/manifest.json)
consolidate eight methods across OrthoBench, corrected QfO and Three Kingdoms.
The exporter uses three hash-pinned retained score reports and two pinned
OrthoBench audit reports; it does not rescore data.
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

## Supplemental Provenance, Version 2

The [five-comparator prediction readback](RETAINED_OB_COMPARATOR_READBACK_20260926.md)
and [upstream evaluator cross-check](RETAINED_OB_UPSTREAM_CROSSCHECK_20260926.md)
are now included as supplemental evidence, not substituted for the original
source rows. Each of those five methods has a prediction record (path, bytes,
SHA-256), parser identities, local and upstream precision/recall/F1, exact RefOG
count, and an explicit false historical-consumption flag in the new manifest.
The original three methods retain their original prediction-provenance objects.

The exporter validates report scopes, exact five-method inventories, content-
identity links to the historical comparison and readback, prediction membership
in both audits, retained metadata, exact counts, and all three metric comparisons.
Hash-and-size links permit report relocation without pretending that original
absolute prediction paths have been relocated or freshly checked. It does not
read raw prediction files or repeat the upstream evaluator. Full native
conversion, command/version/resource provenance remains incomplete.

All 44 focused export/readback/upstream/partition-comparison tests pass;
25 cover the exporter, including inconsistent evidence, partial-attachment
prevention, changed digests, relocation, and output preservation. The new TSV
is byte-identical to the [original TSV](current_benchmark_scores_20260926/scores.tsv).
The entire original export remains unchanged. The OrthoMCL supplement refers
specifically to July predictions: [April/July partition comparison](OB_ORTHOMCL_PARTITION_COMPARISON_20260926.md)
shows differences outside RefOGs despite equal benchmark statistics. Neither
score agreement nor supplemental evidence proves those partitions identical.

```bash
/usr/bin/python3 -S -m benchmark_tools.export_current_benchmark_scores \
  --results benchmark_tools/results --output /tmp/current-benchmark-scores-new
python -m pytest -q tests/unit/test_export_current_benchmark_scores.py
```

No scientific defaults, predictions, uncertainty, timing or native admissions
changed. Full command/version/resource and transitive-provenance consolidation,
controlled timing and the remaining publication requirements are unfinished.
