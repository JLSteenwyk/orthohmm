### Retrospective Native-Cell Comparator Protocol

A separately frozen retrospective analysis compares all eight planned native
cells with full OrthoFinder 3.1.5 and its sequence-only MCL checkpoint, for
SwissTrees F1, precision and recall. Four cells are admitted, leaving 24 of the
48 endpoints estimable; missing endpoints remain in the correction. The 18
common families contain 563 disjoint represented proteins and 10,765 reference
relations. Raw counts use the official halving/unit-prior conversion, equivalent
to PPV=(TP+2)/(TP+FP+4) and TPR=(TP+2)/(TP+FN+4). Every replicate recomputes
harmonic F1 from macro-family precision and recall, not mean family F1.

The calculation uses 100,000 shared multinomial family draws, PCG64 seed
20260920, nominal percentile quantiles 0.025/0.975 and separate 48-endpoint
Bonferroni quantiles 0.05/96 and 1-0.05/96 with linear interpolation. This is
new computation, not interval transplantation from the earlier 24-endpoint
comparator or 42-endpoint internal factorial analysis. Both remain unchanged.
Independent literal-count readback checked statistics, all 24 intervals,
family outcomes and all 48 output rows; deterministic replay is verification,
not another dataset. Development exposure, family exchangeability, merged-
prediction dependence and finite tail resolution limit this conditional
analysis; it is not selection-adjusted or independent confirmation.
[Frozen protocol](NATIVE_QFO_COMPARATOR_UNCERTAINTY_PROTOCOL_20261009.md),
[actual execution/readback](native_qfo_comparator_execution_20261009_v1.json).

