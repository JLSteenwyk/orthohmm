"""Descriptive Pareto comparisons on complete, explicitly named endpoints."""

import math

ENDPOINTS = (
    "OrthoBench: OrthoBench F1", "QfO: VGNC F", "QfO: SwissTrees F",
    "QfO: TreeFam-A F", "QfO: EC", "QfO: GO", "QfO: FAS",
    "ThreeKingdoms: Three Kingdoms F1",
)
PRIMARY_ENDPOINTS = ENDPOINTS[:4]


def metric_frontier(vectors, endpoints=ENDPOINTS):
    required = set(endpoints)
    if not required or len(required) != len(endpoints):
        raise ValueError("Require nonempty unique endpoints")
    complete, excluded = {}, {}
    for label, scores in vectors.items():
        if set(scores) - required:
            raise ValueError("Unknown metric endpoint")
        if any(isinstance(value, bool) or not math.isfinite(value) or not 0 <= value <= 1
               for value in scores.values()):
            raise ValueError("Invalid metric value")
        missing = sorted(required - scores.keys())
        if missing:
            excluded[label] = missing
        else:
            complete[label] = scores
    dominated_by = {label: [] for label in complete}
    for label, scores in complete.items():
        for other, competitor in complete.items():
            if other != label and all(competitor[key] >= scores[key] for key in endpoints) and any(
                    competitor[key] > scores[key] for key in endpoints):
                dominated_by[label].append(other)
    return {"endpoints": list(endpoints),
            "nondominated_complete_labels": sorted(label for label, better in dominated_by.items() if not better),
            "dominated_by": {label: sorted(better) for label, better in dominated_by.items()},
            "excluded_incomplete": excluded,
            "scope": "Development point-estimate dominance only; not uncertainty-aware superiority or independent validation."}
