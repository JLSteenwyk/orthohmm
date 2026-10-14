"""Prepare immutable follow-up inputs after independently admitted pilot results."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json

PROTOCOL = "hillclimb_singleton_tie_followup_20261014.json"
PILOT_REVISION = "1225e3f9b504ee68e90632bf4efdc92e8232703a"


def verified_cell(summary, key, dataset):
    methods = [method for method in summary["methods"] if method["key"] == key]
    if len(methods) != 1:
        raise ValueError("Require a unique admitted configuration")
    cells = [cell for cell in methods[0]["cells"] if cell["dataset"] == dataset
             and cell["status"] == "independent_content_readback_verified"]
    expected = {metric for metric in ENDPOINTS if metric.startswith(dataset + ":")}
    if len(cells) != 1 or set(cells[0]["scores"]) != expected:
        raise ValueError("Require a unique independently admitted complete dataset cell")
    if any(methods[0]["scores"].get(key) != value for key, value in cells[0]["scores"].items()):
        raise ValueError("Admitted method and cell scores disagree")
    return cells[0]


def pilot_gate(summary):
    cell = verified_cell(summary, "singleton_tie:defer_ambiguous", "OrthoBench")
    control = cell["matched_control"]
    if (cell["singleton_tie_policy"] != "defer_ambiguous" or control["singleton_tie_policy"] != "first_hit"
            or cell["source_revision"] != PILOT_REVISION or control["source_revision"] != PILOT_REVISION
            or control["status"] != "independent_content_readback_verified"
            or cell["coverage"]["genes"] != 251378 or cell["coverage"]["species"] != 12
            or control["coverage"]["genes"] != 251378 or control["coverage"]["species"] != 12):
        raise ValueError("Pilot policy, source, control, or complete gene universe differs")
    metric = "OrthoBench: OrthoBench F1"
    a, b = control["scores"][metric], cell["scores"][metric]
    if not 0 <= a <= 1 or not 0 <= b <= 1:
        raise ValueError("Pilot F1 must be finite and between zero and one")
    delta = b - a
    return {"control_f1": a, "variant_f1": b, "delta": delta,
            "eligible": delta >= -.005, "scope": "Prespecified development screen, not equivalence or improvement."}


def prepare(repo, readback, output):
    if output.exists():
        raise FileExistsError("Preserve earlier follow-up preparations")
    summary = json.loads(readback.read_text())
    check(summary["score_table"])
    gate = pilot_gate(summary)
    if not gate["eligible"]:
        raise ValueError("Pilot failed the prespecified development decline screen")
    accounting = subprocess.check_output(["sacct", "-j", "24250", "-X", "-n", "-P",
                                         "-o", "JobID,State%30,ExitCode"], text=True)
    if accounting.strip() != "24250|COMPLETED|0:0":
        raise ValueError("Require successful terminal pilot before follow-up preparation")
    pilot_root = repo / "benchmarks/results/singleton_tie_orthobench_r1_v1_20261014"
    pilot_path = pilot_root / "execution_plan.json"
    pilot_report_path = pilot_root / "results.json"
    pilot = json.loads(pilot_path.read_text())
    pilot_report = json.loads(pilot_report_path.read_text())
    if (pilot["source_revision"] != PILOT_REVISION or pilot_report["job_id"] != "24250"
            or pilot_report["status"] != "both_singleton_policies_scored"
            or pilot_report["execution_plan"] != record(pilot_path)):
        raise ValueError("Admitted pilot differs from the immutable producer")
    rows = pilot_report["arms"]
    if ([row["label"] for row in rows] != ["first_hit", "defer_ambiguous"]
            or any(row["status"] != "native_and_score_independently_verified" for row in rows)
            or rows[0]["scores"]["OrthoBench F1"] != gate["control_f1"]
            or rows[1]["scores"]["OrthoBench F1"] != gate["variant_f1"]):
        raise ValueError("Pilot producer and independent readback disagree")
    protocol = Path(__file__).parent / PROTOCOL
    canonical_root = repo / "benchmarks/results/canonical_downstream_r1_v1_20261014"
    canonical_path = canonical_root / "execution_plan.json"
    canonical = json.loads(canonical_path.read_text())
    original = json.loads((canonical_root / "results.json").read_text())
    science = Path(pilot["source_root"])
    pins = [record(readback), summary["score_table"], record(protocol), record(pilot_path),
            record(pilot_report_path), record(canonical_path), record(canonical_root / "results.json"),
            *summary["reporter_files"], *pilot["sources"], *pilot["runtime"]["files"]]
    cases = []
    for dataset, genes, species in (("QfO", 984137, 78), ("ThreeKingdoms", 443217, 12)):
        control = verified_cell(summary, "canonical_collection:canonical", dataset)
        if (control["coverage"]["genes"], control["coverage"]["species"]) != (genes, species):
            raise ValueError("Follow-up control gene universe differs")
        rows = [row for row in original["arms"] if (row["dataset"], row["label"]) == (dataset, "canonical")]
        configs = [case for case in canonical["cases"] if case["dataset"] == dataset]
        if len(rows) != 1 or len(configs) != 1:
            raise ValueError("Ambiguous original canonical case")
        row, config = rows[0], configs[0]
        checkpoint = Path(row["ordering"]["checkpoint"])
        manifest = row["ordering"]["manifest"]
        if manifest != record(checkpoint / "manifest.json"):
            raise ValueError("Original canonical checkpoint changed")
        native = canonical_root / dataset / "canonical/native"
        golden = {name: record(native / "orthohmm_phylogeny" / name)
                  for name in (*PRODUCTS, "species_tree.rooted.nwk")}
        for name in PRODUCTS:
            expected, actual = control["products"][name], golden[name]
            if (expected["bytes"], expected["sha256"]) != (actual["bytes"], actual["sha256"]):
                raise ValueError("Follow-up native control differs from its admission")
        pins.extend([manifest, *golden.values(), *config["fastas"], *config["scorers"],
                     row["candidate"]["candidate_partition"], row["candidate"]["membership_constraints"],
                     *control["checked"]])
        cases.append({"dataset": dataset, "case": {**config, "baseline_output": str(native)},
                      "checkpoint": str(checkpoint), "checkpoint_sha": manifest["sha256"],
                      "control": control, "original_row_status": row["status"],
                      "candidate": row["candidate"], "golden": golden})
    pins = list({pin["path"]: pin for pin in pins}.values())
    for pin in pins:
        check(pin)
    result = {"status": "prespecified_followup_inputs_prepared_not_submitted", "repo": str(repo),
              "output": str(output), "preparer_source": record(__file__),
              "protocol": record(protocol), "readback": record(readback), "pilot_gate": gate,
              "pilot_accounting": accounting, "scientific_source": str(science),
              "scientific_revision": PILOT_REVISION, "runtime": pilot["runtime"],
              "tool_path": pilot["tool_path"], "cpu": 16, "cases": cases, "inputs": pins,
              "defaults_promoted": False, "goal_completed": False,
              "scope": "Input preparation only. No inference or scoring submitted; execute the fixed protocol with control gates and frozen scientific modules."}
    output.mkdir(parents=True)
    write_json(output / "prepared_inputs.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "readback", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.repo.resolve(), args.readback.resolve(), args.output.resolve())
    print(result["status"])
