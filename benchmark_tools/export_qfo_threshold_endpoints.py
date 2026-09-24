"""Export all native QfO endpoints for the admitted threshold neighborhood."""

import argparse
import csv
import json
import math
from pathlib import Path

from benchmark_tools.audit_qfo_parameter_swiss import validate_admission
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import AXES, validate_records

RESULT_SHA = "0f5dafd4da19a4ed412f33af3ca9122dde37750f82f960712534f917192e292f"
ENDPOINTS = ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")


def scores(assessment):
    native = validate_records(assessment["native_assessments"], assessment["participant"],
                              set(assessment["swiss_reference_families"]))
    if set(assessment["endpoints"]) != set(ENDPOINTS):
        raise ValueError("Incomplete endpoint inventory")
    result = {}
    for endpoint in ENDPOINTS:
        row = assessment["endpoints"][endpoint]
        x, y = [native[endpoint, axis]["metrics"]["value"] for axis in AXES[endpoint]]
        observed = row["native_participant"]
        if (observed["participant_id"] != assessment["participant"]
                or observed["metric_x"] != x or observed["metric_y"] != y
                or (row["axes"]["x_axis"], row["axes"]["y_axis"]) != AXES[endpoint]):
            raise ValueError("Endpoint differs from native assessment")
        value = (2*x*y/(x+y) if x+y else 0.) if endpoint in {"VGNC", "SwissTrees", "TreeFam-A"} else y
        if not math.isclose(row["score"], value, rel_tol=0, abs_tol=1e-12):
            raise ValueError("Endpoint score arithmetic differs")
        result[endpoint] = value
    result["secondary_mean"] = sum(result.values()) / 6
    if not math.isclose(result["secondary_mean"], assessment["secondary_six_metric_mean"], rel_tol=0, abs_tol=1e-12):
        raise ValueError("Secondary mean differs")
    return result


def export(root, output):
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    path = results / "qfo_parameter_uncertainty_threshold_20260923.json"
    uncertainty = read_frozen(path, RESULT_SHA)
    checked = [record(__file__), record(path), *uncertainty["checked_inputs"]]
    for item in checked:
        check(item)
    inv = uncertainty["admission_inventory"]
    inventory = read_frozen(Path(inv["path"]), inv["sha256"])
    rows = []
    for entry in inventory["arms"]:
        row = dict(arm=entry["arm"], status=entry["status"])
        if entry["status"] == "not_admitted":
            row.update({k: None for k in (*ENDPOINTS, "secondary_mean")})
            row["reason"] = entry["reason"]
        else:
            item = entry["admission"]
            report = read_frozen(Path(item["path"]), item["sha256"])
            pair = report["pairs_manifest"]
            conversion = read_frozen(Path(pair["path"]), pair["sha256"])
            validate_admission(entry["arm"], report, conversion)
            row.update(scores(report["assessment"]))
            if not math.isclose(row["SwissTrees"], uncertainty["point_estimates"][entry["arm"]]["F1"], rel_tol=0, abs_tol=1e-6):
                raise ValueError("SwissTrees differs from reconstructed family counts")
            row["reason"] = ""
        rows.append(row)
    for item in checked:
        check(item)
    output.mkdir()
    columns = ["arm", "status", *ENDPOINTS, "secondary_mean", "reason"]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lines = ["# Corrected QfO Threshold Neighborhood", "",
        "| Arm | Status | GO similarity | EC similarity | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | FAS | Secondary mean |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        values = [row[k] for k in (*ENDPOINTS, "secondary_mean")]
        lines.append("| " + " | ".join([row["arm"], row["status"], *["NA" if v is None else f"{v:.6f}" for v in values]]) + " |")
    limitations = ["GO/EC similarity and FAS are not F1; the six-score mean is a project-defined secondary summary.",
        "Missing CPM arms are not zero. These development-exposed results do not promote defaults.",
        "SwissTrees uses rounded native aggregate values here; raw-family uncertainty is retained separately.",
        "No paired uncertainty for the other five endpoints or the secondary mean is established by this table.",
        "FAS sampling is unseeded in the retained scorer; small differences cannot be assigned to parameter changes alone."]
    (output / "scores.md").write_text("\n".join(lines) + "\n\n" + "\n\n".join(limitations) + "\n")
    manifest = dict(status="qfo_threshold_endpoints_exported", rows=rows, checked_records=checked,
        outputs=[record(output / name) for name in ("scores.tsv", "scores.md")],
        limitations=limitations, publication_ready=False)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    export(args.root.resolve(), args.output.absolute())
