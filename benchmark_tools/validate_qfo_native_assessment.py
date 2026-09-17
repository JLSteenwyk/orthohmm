"""Validate native QfO metric identities and aggregation without conflating their semantics."""

import json
import math
from pathlib import Path
import re

from benchmark_tools.qfo_summarize_scores import harmonic_mean

AXES = {"VGNC": ("TPR", "PPV"), "SwissTrees": ("TPR", "PPV"), "TreeFam-A": ("TPR", "PPV"),
        "EC": ("NR_ORTHOLOGS", "avg Schlicker"), "GO": ("NR_ORTHOLOGS", "avg Schlicker"),
        "FAS": ("NR_ORTHOLOGS", "FAS")}


def swiss_families(text):
    # Extract only declaration identities, not Darwin tree expressions or biological truth.
    matches = re.findall(r"^ReconciledTrees\['([^']+)'\] := RecTreeCase\('([^']+)'", text, re.MULTILINE)
    declarations = [line for line in text.splitlines() if line.startswith("ReconciledTrees[")]
    if (not matches or len(matches) != len(declarations) or any(a != b for a, b in matches)
            or len({a for a, b in matches}) != len(matches)):
        raise ValueError("Ambiguous SwissTrees reference declaration inventory")
    return {a for a, b in matches}


def number(value, count=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError("Invalid native QfO numeric value")
    if count and value != int(value):
        raise ValueError("Non-integer assessed relation count")
    return value


def validate_records(records, participant, families):
    expected = {(challenge, metric) for challenge, axes in AXES.items() for metric in axes}
    expected |= {(f"SwissTrees-{family}", metric) for family in families for metric in ("TPR", "PPV")}
    observed, identifiers = {}, set()
    for row in records:
        if (row["type"] != "assessment" or row["community_id"] != "QfO" or row["participant_id"] != participant
                or not isinstance(row["_id"], str) or not row["_id"] or row["_id"] in identifiers):
            raise ValueError("Wrong native assessment identity or duplicate ID")
        identifiers.add(row["_id"])
        key = row["challenge_id"], row["metrics"]["metric_id"]
        if key not in expected or key in observed:
            raise ValueError("Unknown or duplicate native metric")
        value = number(row["metrics"]["value"], count=key[1] == "NR_ORTHOLOGS")
        number(row["metrics"]["stderr"])
        if key[1] != "NR_ORTHOLOGS" and value > 1:
            raise ValueError("QfO score outside unit interval")
        observed[key] = row
    if set(observed) != expected:
        raise ValueError("Missing native challenge or SwissTrees-family metric")
    return observed


def validate_aggregation(data, challenge, participant, native, template):
    if data["type"] != "aggregation" or data["challenge_ids"] != [challenge]:
        raise ValueError("Wrong aggregated challenge identity")
    inline = data["datalink"]["inline_data"]
    axes = inline["visualization"]
    expected_axes = template["datalink"]["inline_data"]["visualization"]
    if axes != expected_axes or (axes["x_axis"], axes["y_axis"]) != AXES[challenge]:
        raise ValueError("Metric axes differ from frozen challenge definition")
    rows = inline["challenge_participants"]
    if len(rows) != 1 or rows[0]["participant_id"] != participant:
        raise ValueError("Unexpected aggregated participant")
    row = rows[0]
    for axis, metric in zip(("x", "y"), AXES[challenge]):
        expected = native[challenge, metric]["metrics"]
        for field, target in (("metric_", "value"), ("stderr_", "stderr")):
            value = number(row[field + axis])
            if value != expected[target]:
                raise ValueError("Aggregation disagrees with native metric or standard error")
    score = harmonic_mean(row["metric_x"], row["metric_y"]) if challenge in {"VGNC", "SwissTrees", "TreeFam-A"} else row["metric_y"]
    return {"score": score, "score_semantics": "harmonic mean of native TPR and PPV" if challenge in {"VGNC", "SwissTrees", "TreeFam-A"}
            else AXES[challenge][1], "native_participant": row, "axes": axes}


def validate_directory(directory, participant, reference):
    reference = Path(reference)
    families = swiss_families((reference / "2020/ReconciledTrees_SwissTrees.drw").read_text())
    native_path = directory / "assessment_out/Assessment_datasets.json"
    records = json.loads(native_path.read_text())
    native = validate_records(records, participant, families)
    endpoints = {}
    paths = [native_path, reference / "2020/ReconciledTrees_SwissTrees.drw"]
    for challenge in AXES:
        path = directory / "results" / challenge / (challenge + ".json")
        template = reference / "data" / (challenge + ".json")
        endpoints[challenge] = validate_aggregation(json.loads(path.read_text()), challenge, participant, native,
                                                    json.loads(template.read_text()))
        paths.extend([path, template])
    return {"participant": participant, "endpoints": endpoints, "native_assessments": records,
            "swiss_reference_families": sorted(families),
            "secondary_six_metric_mean": sum(e["score"] for e in endpoints.values()) / 6,
            "limitations": ["Secondary mean is project-defined, not official QfO F1.",
                            "Native standard errors are not paired method-difference confidence intervals.",
                            "NR_ORTHOLOGS measures challenge-assessed relations, not total submitted pair coverage.",
                            "This validates endpoint content, not scheduler, source/runtime or input provenance."]}, paths
