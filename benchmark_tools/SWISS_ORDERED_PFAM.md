# Ordered-Pfam Annotation Error Analysis

This bounded postprocessing analysis implements the prospective
[protocol](results/SWISS_ORDERED_PFAM_PROTOCOL_20261007.md). It preserves every
corrected SwissTrees family and member. It is not full biological architecture
validation, an independent accuracy evaluation or a new benchmark campaign.

Run each stage once in a fresh namespace, using the exact commit containing
both tested sources. Stage 1 reads retained annotations and alignment lengths,
not prediction counts. Stage 2 independently reconstructs all descriptors with
all-pair interval-overlap checks. Only then may stages 3 and 4 use the existing
54 integer count rows. The score reader uses rational arithmetic, not pooled
counts or mean family F1; it checks both TSVs and every human table row.
JSON/Biopython/csv parsing is shared; numerical and descriptor logic is separate.

```bash
python benchmark_tools/prepare_swiss_ordered_pfam.py features \
  --repo "$REPO" --source-commit "$COMMIT" --output "$FEATURES"
python benchmark_tools/readback_swiss_ordered_pfam.py features \
  --repo "$REPO" --source-commit "$COMMIT" --report "$FEATURES" \
  --output "$FEATURE_READER"
python benchmark_tools/prepare_swiss_ordered_pfam.py scores \
  --repo "$REPO" --source-commit "$COMMIT" --features "$FEATURES" \
  --feature-reader "$FEATURE_READER" --output "$PROJECTION_DIR"
python benchmark_tools/readback_swiss_ordered_pfam.py scores \
  --repo "$REPO" --source-commit "$COMMIT" \
  --report "$PROJECTION_DIR/report.json" --output "$PROJECTION_READER"
```

Use existing Python 3.10.13/Biopython 1.87 with isolated execution (`-I -B`)
and one-thread numerical-library limits for selected invocations. Retain command,
exit status, stderr, source identity and receipts. Any invalid or missing member,
annotation length/coordinate mismatch or incomplete readback fails the complete
stage: no subset scoring or automatic retry. Empty bins remain NA. No annotation
extraction, alignment, tree inference, raw scoring or bootstrap is repeated.

All families are development-exposed. Predicted zero-hit/overlap states do not
establish biological domain absence or complete architecture. Conditional C/R
contrasts are not interactions; initial HMM search remains on with downstream
profiles off. Failed R1 timing remains ineligible; none of these stages provides
new timing, independent confirmation, default promotion or publication readiness.
