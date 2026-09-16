# OrthoBench Uncertainty Protocol

Recorded before running the new bootstrap, after historical aggregate results
were known. This is a development-exposed analysis, not preregistration or
independent confirmation of method selection.

- Retain the complete 70-RefOG benchmark and its low-certainty exclusions.
- Compare retained OrthoHMM high-sensitivity and satellite_v2 phylogenetic
  groups against OrthoFinder 3.1.5 full. No choice between outputs by bootstrap
  outcome. Historical labels do not change this evidence into a held-out test.
- Resample RefOGs jointly across methods, with replacement, 20,000 times;
  NumPy PCG64 seed 20260916. Recompute weighted TP, FP, FN (each divided by
  reference size minus one), then F1, precision, and recall in every draw.
- Report point differences in percentage points and 95% paired percentile
  intervals. Also report Bonferroni-adjusted percentile intervals across all
  six emitted comparisons (two methods times three metrics). These are
  approximate intervals, not exact simultaneous-coverage guarantees.
- Report individual-family F1 wins/ties/losses descriptively. Do not average
  these family F1 values to replace the actual benchmark statistic.
- No pair-level bootstrap, tuning based on intervals, or bootstrap p-values.
- Family exchangeability is an assumption: shared evolutionary histories and
  errors spanning multiple families may create residual dependence. Method
  selection bias remains and requires independent confirmation.

## Output Semantics

The installed OrthoFinder 3.1.5 `file_updates/ogs.py:post_hogs_processing`
reads root hierarchical groups, adds unassigned singletons, and calls
`WriteOrthogroupFiles`. Its `comparative_genomics/orthologues.py` can remove
the redundant root N0 file after finalization. Therefore the retained
`Results_Jul15/Orthogroups/Orthogroups.txt` is not automatically a legacy
sequence-only output merely because of its filename. The run log confirms
version 3.1.5 and successful completion. OrthoHMM uses final root HOGs for
satellite_v2 and `orthogroups_profiles_refined.txt` for high sensitivity.

The analysis JSON records exact input/scorer hashes and per-family sufficient
statistics. Broader scoring and provenance audits remain required for the
publication baseline.
