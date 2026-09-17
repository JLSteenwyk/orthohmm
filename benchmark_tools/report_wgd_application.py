"""Assemble complete biological application diagnostics after native admission."""

from benchmark_tools.bootstrap_wgd_application import intervals
from benchmark_tools.score_wgd_application import score_pairs, summarize

METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "sonicparanoid")
DIAGNOSTIC = "orthofinder_mcl_checkpoint"
CLASSES = ("High", "Low", "Sparse")


def assemble(prepared, references, owners, admitted):
    if set(admitted) != set(METHODS) | {DIAGNOSTIC}:
        raise ValueError("Require every prespecified method and diagnostic status")
    cohort = prepared["cohort_pairs"]
    if any(row["experimental_class"] not in CLASSES for row in cohort):
        raise ValueError("Unrecognized experimental class")
    methods, rows_by_method = {}, {}
    for name in (*METHODS, DIAGNOSTIC):
        entry = admitted[name]
        status = entry["status"]
        if status not in {"admitted", "execution_failed", "admission_failed"}:
            raise ValueError("Unknown or unfinished method status")
        if status != "admitted":
            if entry.get("groups") is not None or not entry.get("reason"):
                raise ValueError("Failed method requires a reason and no fabricated groups")
            methods[name] = {"status": status, "reason": entry["reason"], "summary": None, "rows": None}
            rows_by_method[name] = None
            continue
        rows = score_pairs(cohort, entry["groups"], references, owners)
        rows_by_method[name] = rows
        methods[name] = {"status": "evaluated", "summary": summarize(rows), "rows": rows,
                         "strata": {label: summarize([r for r in rows if r["experimental_class"] == label])
                                    for label in CLASSES}}
    examples = []
    lookup = {name: {tuple(r["orf_pair"]): r for r in rows} if rows is not None else None
              for name, rows in rows_by_method.items()}
    cohort_keys = {tuple(r["orf_pair"]) for r in cohort}
    for example in prepared["prespecified_examples"]:
        key = tuple(example["orf_pair"])
        if key not in cohort_keys:
            raise ValueError("Prespecified example absent from full cohort")
        examples.append({**example, "methods": {name: values[key] if values is not None else None
                                               for name, values in lookup.items()}})
    return {"status": "biological_application_endpoints_assembled", "publication_ready": False,
            "cohort_pairs": len(cohort), "methods": methods, "prespecified_examples": examples,
            "uncertainty": intervals(cohort, {name: rows_by_method[name] for name in METHODS}),
            "limitations": ["Application uses development-exposed Saccharomyces data, not independent generalization.",
                            "Homolog-supported paralog separation is not proof of cross-species copy-specific orthology.",
                            "Foreign-pillar and unmapped members are separate descriptive diagnostics.",
                            "MCL checkpoint is diagnostic, not a separately run sequence-only OrthoFinder pipeline.",
                            "Caller must verify native admission and artifact identities before using this assembly."]}
