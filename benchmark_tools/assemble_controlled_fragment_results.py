"""Admit terminal fragment outputs and report every fixed identity and failure."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_simulation_results import TERMINAL, admit_method, input_universe
from benchmark_tools.controlled_fragment_observations import CONDITION, SEEDS, METHODS, METRICS, paired_intervals, score_strata
from benchmark_tools.prepare_controlled_fragment_observations import PINS, PROTOCOL, record
from benchmark_tools.run_controlled_fragment_methods import observed_inputs
from benchmark_tools.run_simulation_generation import verify_file
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment
from benchmark_tools.simulation_method_outputs import load_predictions


def terminal_job(accounting, job):
    rows = [row for row in csv.DictReader(io.StringIO(accounting), delimiter="|") if row.get("JobIDRaw") == str(job)]
    if len(rows) != 1 or rows[0]["State"].split()[0] not in TERMINAL:
        raise ValueError("Fragment job must be uniquely terminal; do not score a partial run")
    if rows[0]["State"] == "COMPLETED" and rows[0]["ExitCode"] != "0:0":
        raise ValueError("Contradictory scheduler completion")
    return rows[0]


def prediction_scores(method, dataset, execution, runtime, verified, owners, species, truth, flags):
    admission = admit_method(method, dataset, execution, verified, runtime)
    if admission["status"] != "admitted":
        return {**admission, "status": "failed"}
    pairs, paths = load_predictions(method, Path(dataset["methods"][method]["output"]), owners, species)
    parent = "orthofinder_full" if method == "orthofinder_sequence_only" else method
    inventory = {Path(item["absolute_path"]).resolve() for item in execution["methods"][parent]["outputs"]}
    if not {p.resolve() for p in paths} <= inventory:
        raise ValueError("Fragment predictions not in native execution inventory")
    return {"status": "complete", "admission": admission, "prediction_artifacts": [record(p) for p in paths],
            **score_strata(pairs, truth["ortholog_pairs"], owners, flags)}


def fragment_rows(dataset, panel_row, provenance, runtime):
    common = {"arm": "fragment", "seed": dataset["seed"], "condition": CONDITION}
    if panel_row is None or "execution" not in panel_row:
        reason = panel_row.get("reason", panel_row["status"]) if panel_row else "Dataset was not reached before terminal interruption"
        return [{**common, "method": method, "status": "unavailable", "failure_stage": "execution_unavailable", "reason": reason} for method in METHODS]
    item = panel_row["execution"]
    status_path = Path(item["absolute_path"])
    expected_status = Path(dataset["methods"][METHODS[0]]["output"]).parent / "execution/status.json"
    if status_path.resolve() != expected_status.resolve():
        raise ValueError("Fragment status belongs to a different inference identity")
    verify_file(status_path, item)
    execution = json.loads(status_path.read_text())
    verified = observed_inputs(dataset)
    if (execution["dataset"] != dataset["label"] or execution["verified_inputs"] != verified
            or execution["provenance"] != provenance or execution["status"] != "finished_pending_native_validation"):
        raise ValueError("Fragment execution identity, terminal state or provenance differs")
    truth = json.loads(Path(dataset["truth"]).read_text())
    owners, species = input_universe(verified["inputs"], truth)
    coordinates = json.loads(Path(dataset["coordinates"]["absolute_path"]).read_text())
    flags = {row["gene"]: row["fragment"] for row in coordinates}
    rows = []
    for method in METHODS:
        try:
            outcome = prediction_scores(method, dataset, execution, runtime, verified, owners, species, truth, flags)
        except ValueError as error:
            outcome = {"status": "unavailable", "failure_stage": "native_integrity_or_conversion", "reason": str(error)}
        parent = "orthofinder_full" if method == "orthofinder_sequence_only" else method
        evidence = status_path.parent
        resource_paths = [p for p in (evidence / (parent + ".time.log"), evidence / (parent + ".log")) if p.is_file()]
        rows.append({**common, "method": method, "execution": item, "truth_sha256": verified["truth"]["sha256"],
                     "resource_evidence": [record(p) for p in resource_paths],
                     "independent_timing": method != "orthofinder_sequence_only", **outcome})
    return rows


def metric_means(records):
    summaries = []
    for arm in ("baseline", "fragment"):
        for method in METHODS:
            rows = [r for r in records if r["arm"] == arm and r["method"] == method]
            if len(rows) != len(SEEDS) or {r["seed"] for r in rows} != set(SEEDS):
                raise ValueError("Incomplete seed inventory for metric summaries")
            row = {"arm": arm, "method": method, "planned_seeds": len(SEEDS),
                   "successful_seeds": [r["seed"] for r in rows if r["status"] == "complete"],
                   "failed_or_unavailable_seeds": [r["seed"] for r in rows if r["status"] != "complete"], "metrics": {}}
            for metric in METRICS:
                eligible = [r for r in rows if r["status"] == "complete" and r["score"][metric] is not None]
                row["metrics"][metric] = {"mean": sum(r["score"][metric] for r in eligible) / len(eligible) if eligible else None,
                                          "eligible_seeds": [r["seed"] for r in eligible]}
            summaries.append(row)
    return summaries


def tables(output, records, comparisons):
    fields = ["arm", "seed", "method", "status", "fragment_endpoints", "tp", "fp", "fn", "f1", "precision", "recall", "predicted_pairs", "eligible_true_pairs", "pair_endpoint_coverage", "reason"]
    for name, stratified in (("scores.tsv", False), ("strata.tsv", True)):
        with (output / name).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in records:
                strata = row.get("strata", []) if stratified else [{"fragment_endpoints": "all", "score": row.get("score", {})}]
                if stratified and not strata:
                    strata = [{"fragment_endpoints": count, "score": {}} for count in range(3)]
                for stratum in strata:
                    value = {k: row[k] for k in ("arm", "seed", "method", "status")}
                    value.update(fragment_endpoints=stratum["fragment_endpoints"], reason=row.get("reason", ""))
                    value.update({key: stratum["score"].get(key) for key in fields[5:-1]})
                    writer.writerow({k: "NA" if v is None else v for k, v in value.items()})
    fields = ["target_arm", "target_method", "reference_arm", "reference_method", "metric", "status", "paired_seeds", "estimate", "nominal_low", "nominal_high", "adjusted_low", "adjusted_high"]
    with (output / "comparisons.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in comparisons:
            value = {k: row[k] for k in fields[:6]}
            value.update(paired_seeds=len(row["eligible_seeds"]), estimate=row["estimate"])
            for kind in ("nominal", "adjusted"):
                interval = row[kind + "_interval"]
                value.update({kind + "_low": interval[0] if interval else None, kind + "_high": interval[1] if interval else None})
            writer.writerow({k: "NA" if v is None else v for k, v in value.items()})


def assemble(root, manifest_path, manifest_sha, runtime_path, runtime_sha, execution_root, job, output):
    if output.exists():
        raise FileExistsError("Existing fragment results; do not overwrite scientific evidence")
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS"], text=True)
    scheduler = terminal_job(accounting, job)
    panel = read_frozen(manifest_path, manifest_sha)
    if (panel.get("schema") != "controlled_fragment_observations_v1" or panel["status"] != "prepared_unexecuted"
            or [d["seed"] for d in panel["datasets"]] != list(SEEDS) or panel["inference_identities"] != 30):
        raise ValueError("Wrong prepared fragment inventory")
    for item in panel["baseline_pins"]:
        if PINS.get(Path(item["absolute_path"]).name) != item["sha256"]:
            raise ValueError("Unknown fragment parent pin")
        verify_file(Path(item["absolute_path"]), item)
    if panel["protocol"] != record(root / "benchmark_tools/results" / PROTOCOL):
        raise ValueError("Fragment protocol changed")
    runtime = read_frozen(runtime_path, runtime_sha)
    verify_environment(runtime)
    status_path = execution_root / "panel_status.json"
    status = json.loads(status_path.read_text())
    provenance = status["provenance"]
    if (status["schema"] != "controlled_fragment_panel_execution_v1"
            or status["status"] not in {"finished_pending_native_validation", "interrupted_or_integrity_failure"}
            or provenance["fragment_manifest"] != record(manifest_path)
            or provenance["runtime_manifest"] != record(runtime_path)
            or provenance["allocation"]["job_id"] != str(job)
            or status["accuracy_evaluated"] is not False):
        raise ValueError("Fragment panel execution identity or terminal status differs")
    for item in provenance["sources"]:
        verify_file(Path(item["absolute_path"]), item)
    execution_rows = {row["label"]: row for row in status["datasets"]}
    if len(execution_rows) != len(status["datasets"]) or set(execution_rows) - {d["label"] for d in panel["datasets"]}:
        raise ValueError("Duplicate or unknown executed fragment dataset")
    report = {"schema": "controlled_fragment_results_v1", "status": "admitting", "publication_ready": False,
              "manifest": record(manifest_path), "runtime": record(runtime_path), "execution": record(status_path),
              "source": record(__file__), "scheduler": scheduler, "accounting_raw": accounting,
              "records": [], "limitations": ["Development-exposed synthetic truncations, not natural fragment truth or independent biology.",
                  "Whole-seed percentile intervals are conditional approximations from at most ten seeds; failure exclusions are explicit.",
                  "Fixed15-endpoint Bonferroni quantiles do not establish exact simultaneous coverage.",
                  "Current OrthoHMM inventory differs from historical controls; required scientific versions and frozen method sources remain unchanged, but historical output equivalence is unproven.",
                  "Shared-host resource observations are descriptive, not isolated speed rankings or paired baseline timings."]}
    output.mkdir(parents=True, exist_ok=False)
    try:
        for dataset in panel["datasets"]:
            observed_inputs(dataset)
            baseline = dataset["baseline_records"]
            if len(baseline) != 4 or {r["method"] for r in baseline} != set(METHODS) or any(r["seed"] != dataset["seed"] or r["arm"] != "baseline" or r["status"] != "complete" for r in baseline):
                raise ValueError("Malformed retained baseline score inventory")
            report["records"].extend(baseline)
            report["records"].extend(fragment_rows(dataset, execution_rows.get(dataset["label"]), provenance, runtime))
        report["comparisons"] = paired_intervals(report["records"])
        report["metric_means"] = metric_means(report["records"])
        tables(output, report["records"], report["comparisons"])
        report["status"] = "complete_with_explicit_outcomes"
        report["fragment_successes"] = sum(r["arm"] == "fragment" and r["status"] == "complete" for r in report["records"])
    except Exception as error:
        report.update(status="assembly_failed", error_type=type(error).__name__, reason=str(error))
        raise
    finally:
        (output / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--runtime-sha256", required=True)
    parser.add_argument("--execution-root", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = assemble(args.root.resolve(), args.manifest.resolve(), args.manifest_sha256,
                      args.runtime.resolve(), args.runtime_sha256, args.execution_root.resolve(), args.job, args.output.resolve())
    print(json.dumps({"status": result["status"], "records": len(result["records"]), "fragment_successes": result["fragment_successes"], "comparisons": len(result["comparisons"])}))
