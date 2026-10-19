"""Independently admit the complete controlled projection/HMMER recovery panel."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Bio import SeqIO

from benchmark_tools import report_profile_evidence_pilot as base
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_profile_weighting_pilot import advancement, experiment_settings
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, validate_ownership
from orthohmm.accuracy import load_accuracy_checkpoint


PROTOCOL = "hmm_projection_recovery_readback_20261019.json"
LEDGER = "benchmark_tools/results/hmm_evidence_campaign_status_20261018.json"
UNIVERSE = (251378, 12)
ARMS = experiment_settings("projection_recovery")[0]


def launch(repo, revision, *, entry_key="projection_recovery_pilot"):
    if not isinstance(revision, str) or re.fullmatch(r"[0-9a-f]{40}", revision) is None:
        raise ValueError("Require an explicit immutable full launch revision")
    if subprocess.check_output(["git", "cat-file", "-t", revision], cwd=repo, text=True).strip() != "commit":
        raise ValueError("Launch revision is not a commit")
    blob = subprocess.check_output(["git", "show", revision + ":" + LEDGER], cwd=repo)
    if entry_key not in ("projection_recovery_pilot", "stockholm_annotation_recovery_pilot"):
        raise ValueError("Unknown explicit recovery launch entry")
    entry = json.loads(blob)[entry_key]
    inference = entry["inference"]
    if (inference.get("submitted") is not True or not str(inference.get("job_id", "")).isdigit()
            or tuple(entry["arms"]) != ARMS or inference.get("planned_cpus") != 16
            or inference.get("planned_memory") != "128G"
            or inference.get("planned_wall_limit") != "1-00:00:00"):
        raise ValueError("Require the submitted original recovery job and unchanged planned panel/resources")
    reader = entry["independent_reader"]
    if (re.fullmatch(r"[0-9a-f]{40}", reader.get("runner_commit", "")) is None
            or Path(reader["scientific_reader_source"]).is_absolute()
            or not (repo / reader["scientific_reader_source"]).resolve().is_relative_to(repo.resolve())):
        raise ValueError("Recovery launch must bind the separately frozen reader revision and repository path")
    pin = {key: entry["execution_plan"][key] for key in ("path", "bytes", "sha256")}
    pin["path"] = str((repo / pin["path"]).resolve())
    binding = {"revision": revision, "path": LEDGER, "git_blob_sha256": hashlib.sha256(blob).hexdigest()}
    if entry_key != "projection_recovery_pilot": binding["entry_key"] = entry_key
    return entry, pin, binding


def point(cell):
    return {"f1": cell["scores"][base.ENDPOINTS[0]], "precision": cell["details"]["precision"],
            "recall": cell["details"]["recall"]}


def signature(pins):
    return {key: {field: pin[field] for field in ("bytes", "sha256")} for key, pin in pins.items()}


def verify_panel(cells, rows):
    if tuple(c["label"] for c in cells) != ARMS or tuple(r["label"] for r in rows) != ARMS:
        raise ValueError("Require the complete ordered five-arm recovery panel")
    for i, cell in enumerate(cells):
        expected_policy = "legacy" if i == 0 else "residue_preserving"
        if (cell.get("status") != "independent_content_readback_verified" or cell.get("dataset") != "OrthoBench"
                or cell["activation"].get("projection_policy") != expected_policy
                or cell["activation"].get("arm") != ARMS[i]):
            raise ValueError("Unadmitted recovery cell or wrong profile projection policy")
    legacy, projection, builtin, measurement, hmmer = cells
    if len({cell["activation"]["candidate_digest"] for cell in cells[2:]}) != 1:
        raise ValueError("Recovery active candidate nominations differ")
    if (measurement["graph_stages"] != projection["graph_stages"]
            or signature(measurement["products"]) != signature(projection["products"])
            or signature({"tree": measurement["rooted_species_tree"]}) != signature({"tree": projection["rooted_species_tree"]})
            or signature({k: measurement["candidate"][k] for k in ("candidate_partition", "membership_constraints")})
            != signature({k: projection["candidate"][k] for k in ("candidate_partition", "membership_constraints")})
            or measurement["scores"] != projection["scores"] or measurement["details"] != projection["details"]
            or measurement["coverage"] != projection["coverage"]
            or rows[3].get("measurement_control_matches") is not True):
        raise ValueError("HMMER measurement differs from the repaired control")
    for i in (1, 2, 4):
        cell, row = cells[i], rows[i]
        gate = advancement(point(legacy), point(cell)) if i == 1 else base.selective_gate(point(legacy), cell)
        if gate != row["advancement"]:
            raise ValueError("Original legacy-relative recovery gate differs")
        cell["advancement"] = gate
        if i != 1:
            matched = advancement(point(projection), point(cell))
            if (matched != row["matched_projection_advancement"]
                    or row.get("accuracy_interpretation") != cell["activation"]["applicability"]):
                raise ValueError("Matched repaired-control gate or applicability differs")
            cell["matched_projection_advancement"] = matched
    if "advancement" in rows[0] or "advancement" in rows[3]:
        raise ValueError("Legacy and measurement controls cannot become advancing candidates")
    return {"f1": point(hmmer)["f1"] - point(builtin)["f1"],
            "precision": point(hmmer)["precision"] - point(builtin)["precision"],
            "recall": point(hmmer)["recall"] - point(builtin)["recall"],
            "scope": "HMMER-minus-builtin under matched repaired construction and decision rules. Legacy-relative gains alone cannot establish HMMER-specific benefit; development point contrast, not independent confirmation."}


def uncertainty(cells):
    spec = json.loads(Path(__file__).with_name(PROTOCOL).read_text())
    expected = {"replicates": 100000, "seed": 20261018, "alpha": .05, "families": 70,
                "metrics": list(base.METRICS), "multiplicity_endpoints": 24}
    if (tuple(spec["ordered_arms"]) != ARMS or tuple(c["label"] for c in cells) != ARMS
            or any(spec["uncertainty"][key] != value for key, value in expected.items())):
        raise ValueError("Changed recovery uncertainty contrasts or multiplicity")
    return {"versus_legacy": base.uncertainty(cells, multiplicity_endpoints=24),
            "versus_projection": base.uncertainty(cells[1:], multiplicity_endpoints=24),
            "hmmer_versus_builtin": base.uncertainty([cells[2], cells[4]], multiplicity_endpoints=24),
            "scope": "All24 prespecified metric contrasts; descriptive paired-family intervals on development-exposed data. No change to point gates, pair-independent inference or independent biological confirmation."}


def ownership(plan, report):
    names, species, q, t, scores = load_accuracy_checkpoint(Path(plan["case"]["checkpoint"]), verify=True)
    del q, t, scores
    owners = {}
    for pin in plan["case"]["fastas"]:
        check(pin)
        for sequence in SeqIO.parse(pin["path"], "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty recovery input sequence")
            owners[sequence.id] = Path(pin["path"]).name
    actual = json.loads(json.dumps(validate_ownership(names, species, owners)))
    if actual != report["ownership"] or (actual["genes"], actual["species"]) != UNIVERSE:
        raise ValueError("Recovery checkpoint/FASTA universe or ownership differs")
    return names


def reporter_inventory(driver):
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    pins = sources(driver) + [record(driver / "benchmark_tools" / name) for name in (
        PROTOCOL, base.PROTOCOL, "hmm_projection_integrity_recovery_20261018.json")]
    for pin in pins:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Recovery reader differs from committed source")
    return revision, pins


def next_actions(cells):
    result = []
    for i in (1, 2, 4):
        cell = cells[i]
        if i != 1 and cell["activation"]["completely_measured_queries"] == 0:
            action = "Original applicability fallback: probe corrected QfO then Three Kingdoms without retuning; incomplete measurements are not negative accuracy."
        elif cell["advancement"]["eligible"]:
            action = "Run unchanged corrected QfO and Three Kingdoms, retain all eight endpoints/populations/costs and complete promising production opt-ins; no default promotion from a pilot."
        else:
            action = "Stop this fixed branch without retuning or full-vector extension; retain negative evidence."
        result.append({"method": cell["label"], "action": action})
    result.append({"method": "campaign", "action": "Retain completed weighting/local/alignment evidence, assess justified combinations only after individual effects, freeze candidate before confirmation or document unavailable independent data, and finish final comparison/recommendation/remote audit."})
    return result


def run(repo, revision, output, *, entry_key="projection_recovery_pilot"):
    if output.exists():
        raise FileExistsError("Preserve every recovery readback attempt")
    entry, plan_pin, binding = (launch(repo, revision) if entry_key == "projection_recovery_pilot"
                                else launch(repo, revision, entry_key=entry_key))
    check(plan_pin)
    plan = json.loads(Path(plan_pin["path"]).read_text())
    root, source = Path(plan["output"]), Path(plan["source_root"])
    if (Path(plan_pin["path"]).parent != root or Path(plan["repo"]) != repo
            or source != (repo / entry["scientific_source"]).resolve()
            or plan["source_revision"] != entry["runner_commit"]):
        raise ValueError("Recovery plan differs from immutable launch")
    producer = record(root / "results.json")
    report = json.loads(Path(producer["path"]).read_text())
    if report["execution_plan"] != plan_pin:
        raise ValueError("Recovery producer uses a different execution plan")
    base.require_completed(report, plan, ARMS, "projection_recovery", entry["inference"]["job_id"])
    accounting = base.terminal(entry["inference"]["job_id"])
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != plan["source_revision"]:
        raise ValueError("Frozen recovery scientific source moved")
    if (sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")] != plan["sources"]
            or runtime() != plan["runtime"]):
        raise ValueError("Recovery scientific inventory or runtime changed")
    checked = [plan_pin, producer, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in checked:
        check(pin)
    names = ownership(plan, report)
    driver = Path(__file__).resolve().parents[1]
    reader_revision, reader_files = reporter_inventory(driver)
    if (reader_revision != entry["independent_reader"]["runner_commit"]
            or driver != (repo / entry["independent_reader"]["scientific_reader_source"]).resolve()):
        raise ValueError("Actual recovery reader differs from the immutable launch binding")
    context, context_pins = base.retained_context(repo)
    checked += [*reader_files, *context_pins]
    output.mkdir(parents=True)
    result = {"status": "running", "experiment": "projection_recovery", "launch_ledger": binding,
        "execution_plan": plan_pin, "producer": producer, "scheduler": accounting,
        "source_revision": plan["source_revision"], "reporter_revision": reader_revision,
        "reporter_files": reader_files, "protocol": record(driver / "benchmark_tools" / PROTOCOL),
        "context": context, "cells": [], "accuracy_result_admitted": False,
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Separate shared-helper source/native/score readback, not independent biological confirmation. Development-exposed pilot; seven endpoint slots remain unrun until gated extensions. Shared-host cached incremental timings and max-process RSS are not fresh isolated full-pipeline costs. Functional means are not F1."}
    write_json(output / "summary.json", result)
    try:
        for row in report["arms"]:
            projection = result["cells"][1] if len(result["cells"]) >= 2 else None
            cell = base.admit_arm(repo, plan, row, "projection_recovery", names, projection_control=projection)
            result["cells"].append(cell)
            checked += cell["checked"]
            write_json(output / "summary.json", result)
        result["backend_contrast"] = verify_panel(result["cells"], report["arms"])
        result["uncertainty"] = uncertainty(result["cells"])
        result["required_next_actions"] = next_actions(result["cells"])
        base.write_table(output / "scores.tsv", result["cells"], context)
        result["score_table"] = record(output / "scores.tsv")
        checked.append(result["score_table"])
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Runtime changed during recovery readback")
        result["checked"] = list({pin["path"]: pin for pin in checked}.values())
        result["advancing_labels"] = [c["label"] for c in result["cells"] if c.get("advancement", {}).get("eligible")]
        result.update(status="projection_recovery_independently_read_back", accuracy_result_admitted=True)
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--launch-revision", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--launch-entry", choices=("projection_recovery_pilot", "stockholm_annotation_recovery_pilot"),
                        default="projection_recovery_pilot")
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.launch_revision, args.output.resolve(), entry_key=args.launch_entry)["status"])
