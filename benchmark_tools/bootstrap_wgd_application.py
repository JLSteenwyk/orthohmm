"""Frozen paired pillar bootstrap for the descriptive WGD application endpoints."""

import numpy as np

ENDPOINTS = ("separation_rate", "supported_separation_rate", "mean_non_scer_coverage")
CONTRASTS = (("orthohmm_satellite_v2", "orthohmm_high_sensitivity"),
             ("orthohmm_satellite_v2", "orthofinder_full"),
             ("orthohmm_satellite_v2", "sonicparanoid"),
             ("orthohmm_high_sensitivity", "orthofinder_full"))


def intervals(cohort, method_rows):
    population = sorted((r for r in cohort if r["reference_eligible"]),
                        key=lambda r: (r["reference_pillar"], tuple(r["orf_pair"])))
    keys = [tuple(r["orf_pair"]) for r in population]
    if not keys or len(set(keys)) != len(keys):
        raise ValueError("Empty or duplicated fixed comparison population")
    pillars = sorted({r["reference_pillar"] for r in population})
    unit = {p: i for i, p in enumerate(pillars)}
    pair_units = np.array([unit[r["reference_pillar"]] for r in population])
    coverage_denominators = np.array([sum(n for s, n in r["available_members_by_species"].items()
                                         if s != "Scerevisiae") for r in population])
    values = {}
    for method in {m for pair in CONTRASTS for m in pair}:
        rows = method_rows.get(method)
        if rows is None:
            continue
        mapping = {tuple(r["orf_pair"]): r for r in rows}
        if len(mapping) != len(rows) or set(mapping) != {tuple(r["orf_pair"]) for r in cohort}:
            raise ValueError("Method rows differ from complete experimental cohort")
        values[method] = {}
        for endpoint in ENDPOINTS:
            column = []
            for index, pair in enumerate(population):
                row = mapping[tuple(pair["orf_pair"])]
                if not row["reference_eligible"] or row["reference_pillar"] != pair["reference_pillar"]:
                    raise ValueError("Method changed fixed reference population")
                value = row[endpoint]
                if row["coverage_denominator"] != int(coverage_denominators[index]):
                    raise ValueError("Method-dependent coverage denominator")
                if endpoint == "mean_non_scer_coverage" and coverage_denominators[index] == 0:
                    if value is not None:
                        raise ValueError("Zero denominator must be unavailable")
                    column.append(np.nan)
                elif value is None or not np.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError("Invalid or outcome-dependent missing endpoint")
                elif endpoint != "mean_non_scer_coverage" and value not in (0, 1):
                    raise ValueError("Separation endpoints must be binary")
                else:
                    column.append(value)
            values[method][endpoint] = np.array(column, dtype=float)

    draws = np.random.Generator(np.random.PCG64(20260920)).integers(0, len(pillars), (20000, len(pillars)))
    comparisons = []
    for first, second in CONTRASTS:
        for endpoint in ENDPOINTS:
            result = {"first": first, "second": second, "endpoint": endpoint}
            if first not in values or second not in values:
                comparisons.append({**result, "status": "unavailable_method", "difference_pp": None})
                continue
            delta = values[first][endpoint] - values[second][endpoint]
            valid = np.isfinite(delta)
            totals = np.bincount(pair_units[valid], weights=delta[valid], minlength=len(pillars))
            counts = np.bincount(pair_units[valid], minlength=len(pillars))
            sample_n = counts[draws].sum(axis=1)
            undefined = sample_n == 0
            bootstrap = np.divide(totals[draws].sum(axis=1), sample_n,
                                  out=np.full(20000, np.nan), where=~undefined) * 100
            result.update(pairs=int(valid.sum()), pillars=len(pillars), undefined_replicates=int(undefined.sum()),
                          difference_pp=float(delta[valid].mean() * 100) if valid.any() else None,
                          wins=int((delta[valid] > 0).sum()), losses=int((delta[valid] < 0).sum()),
                          ties=int((delta[valid] == 0).sum()), nominal95_pp=None, bonferroni12_pp=None,
                          status="undefined_replicates" if undefined.any() else "evaluated")
            if not undefined.any():
                result["nominal95_pp"] = np.quantile(bootstrap, [0.025, 0.975], method="linear").tolist()
                result["bonferroni12_pp"] = np.quantile(bootstrap, [0.05 / 24, 1 - 0.05 / 24], method="linear").tolist()
            comparisons.append(result)
    return {"seed": 20260920, "rng": "PCG64", "replicates": 20000, "multiplicity": 12,
            "population_pairs": len(population), "population_pillars": len(pillars),
            "comparisons": comparisons,
            "interpretation": "Exploratory paired differences on development-exposed data; not orthology F1 or independent superiority."}
