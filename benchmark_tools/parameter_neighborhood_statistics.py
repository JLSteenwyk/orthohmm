"""Fixed six-variant uncertainty summary after native/output admission."""

from benchmark_tools.bootstrap_orthobench import METRICS, paired_bootstrap, statistics, weighted_records

BASELINE = "control"
VARIANTS = ("cpm_low", "cpm_high", "norm_low", "norm_high", "margin_low", "margin_high")


def summarize(records, failures):
    """Consume admitted sufficient statistics and externally verified terminal failures.

    This statistical component does not itself admit native runs. Its caller must
    verify output provenance and official-score equivalence before supplying data.
    """
    if BASELINE not in records or BASELINE in failures:
        raise ValueError("A successfully scored control is required")
    if set(records) & set(failures) or set(records) | set(failures) != {BASELINE, *VARIANTS}:
        raise ValueError("Account for exactly all six planned variants once")
    if any(not isinstance(reason, str) or not reason.strip() for reason in failures.values()):
        raise ValueError("Every terminal failure requires a documented reason")
    if len(records) > 1:
        result = paired_bootstrap(records, BASELINE, replicates=20000, seed=20260918,
                                  multiplicity_endpoints=18)
    else:
        names, _, counts = weighted_records(records[BASELINE])
        result = {"baseline": BASELINE, "families": names, "replicates": 0,
                  "planned_replicates": 20000, "seed": 20260918, "alpha": .05,
                  "point_estimates_percent": {BASELINE: dict(zip(METRICS, statistics(counts.sum(axis=0)).tolist()))},
                  "comparisons": {}, "multiplicity": "Bonferroni tail adjustment over 18 planned endpoints",
                  "limitations": ["All six variants failed; no bootstrap intervals were computed."]}
    result.update(planned_variants=list(VARIANTS), planned_endpoints=18,
                  failed_variants=dict(failures), publication_ready=False)
    result["missing_intervals"] = {name: {metric: {"status": "not_estimable_failed_run", "reason": failures[name]}
                                           for metric in METRICS} for name in VARIANTS if name in failures}
    result["limitations"].extend([
        "Post-development prespecified parameter sensitivity, not independent confirmation or default selection.",
        "All six variants and 18 endpoints remain in the planned family even when an inference run fails.",
        "Failed runs have no accuracy estimate or confidence interval; they are not assigned zero scores.",
        "Intervals including zero do not demonstrate equivalence or general robustness."])
    return result
