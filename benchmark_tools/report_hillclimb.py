"""Read back completed hill-climb cells and report separate benchmark endpoints."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo, SeqIO
from orthohmm.files import fetch_fasta_files
from benchmark_tools.admit_qfo_factorial_cell import check_pairs
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.metric_frontier import ENDPOINTS, PRIMARY_ENDPOINTS, metric_frontier
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import validate_confidence, validate_native
from benchmark_tools.run_membership_hillclimb import (
    DATASETS, PRODUCTS, argument, compare_golden, score_orthobench, write_json,
)
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.summarize_three_kingdoms_parity import parse_score
from benchmark_tools.validate_qfo_native_assessment import validate_directory

METRIC_NAMES = {
    "OrthoBench": {"weighted_refog_F1": "OrthoBench F1"},
    "QfO": {"VGNC": "VGNC F", "SwissTrees": "SwissTrees F", "TreeFam-A": "TreeFam-A F",
            "EC": "EC", "GO": "GO", "FAS": "FAS"},
    "ThreeKingdoms": {"BUSCO_reference_pair_F1": "Three Kingdoms F1"},
}
EXPECTED_GENES = {"OrthoBench": 251378, "QfO": 984137, "ThreeKingdoms": 443217}
ROUND_ARMS = {
    "membership": ("control_direct", "control_path", "broader_direct", "broader_path"),
    "full_search": ("control100", "cap400"),
}


def reporter_sources():
    names = ("metric_frontier.py", "run_full_search_hillclimb.py", "run_membership_hillclimb.py",
             "score_ygob_groups.py", "validate_qfo_native_assessment.py", "admit_qfo_factorial_cell.py",
             "audit_historical_profile_ablation.py", "qfo_summarize_scores.py", "simulation_method_outputs.py",
             "summarize_three_kingdoms_parity.py", "orthobench_stage_diagnostics.py", "score_orthobench_partition.py",
             "report_deferred_cells.py", "run_deferred_copy_split.py", "deferred_qfo_inferred.py",
             "run_deferred_inferred_supplement.py", "run_held_tree_hillclimb.py",
             "run_three_kingdoms_membership_supplement.py", "prepare_ob_candidate_neighborhood.py",
             "run_qfo_recovered_assessment.py")
    return [record(Path(__file__)), *[record(Path(__file__).with_name(name)) for name in names]]


def endpoint_vector(dataset, scores):
    expected = set(METRIC_NAMES[dataset].values())
    if set(scores) != expected:
        raise ValueError("Missing or unknown dataset endpoint")
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not 0 <= v <= 1
           for v in scores.values()):
        raise ValueError("Invalid dataset metric value")
    vector = {dataset + ": " + metric: value for metric, value in scores.items()}
    metric_frontier({"validation": vector})
    return vector


def retained_methods(registry):
    methods, seen = {}, set()
    for row in registry["rows"]:
        dataset, key = row["dataset"], row["key"]
        if (dataset, key) in seen:
            raise ValueError("Duplicate retained method/dataset cell")
        seen.add((dataset, key))
        mapping = METRIC_NAMES[dataset]
        if set(row["scores"]) != set(mapping):
            raise ValueError("Retained endpoint schema differs")
        scores = endpoint_vector(dataset, {mapping[k]: v for k, v in row["scores"].items()})
        method = methods.setdefault(key, {"key": key, "label": row["label"], "scores": {}, "cells": []})
        method["scores"].update(scores)
        method["cells"].append({"dataset": dataset, "status": "retained_registered_result",
            "version": row["declared_version"], "semantics": row["output_semantics"],
            "resources": row["resources"], "limitations": row["gaps"],
            "retained_details": {k: row[k] for k in ("counts", "score_details", "submitted_pairs", "retained_pairs", "retained_ob_metadata") if k in row}})
    return methods


def snapshot(path, directory, name):
    data = Path(path).read_bytes()
    destination = directory / name
    destination.write_bytes(data)
    return json.loads(data), {"snapshot": record(destination), "source_path": str(Path(path).resolve()),
        "source_bytes": len(data), "source_sha256": hashlib.sha256(data).hexdigest()}


def replay_native(row, case):
    outputs = {Path(item["path"]).name: item for item in row["outputs"]}
    if not set(PRODUCTS) <= outputs.keys() or len(outputs) != len(row["outputs"]):
        raise ValueError("Missing or duplicate native products")
    parents = {Path(outputs[name]["path"]).parent for name in PRODUCTS}
    if len(parents) != 1:
        raise ValueError("Mixed native output directories")
    phylo = parents.pop()
    argv = case["command"]
    if Path(argument(argv, "--output-directory")) / "orthohmm_phylogeny" != phylo:
        raise ValueError("Native output differs from planned directory")
    fasta = Path(argument(argv, "--fasta-directory"))
    provenance = json.loads((phylo / "provenance_manifest.json").read_text())
    expected = {"root_duplication_rule": argument(argv, "--root-rule"),
                "pair_orthology_rule": argument(argv, "--pair-rule"),
                "species_tree_mode": argument(argv, "--species-tree-mode"),
                "species_tree_rooting": argument(argv, "--species-tree-rooting")}
    if any(provenance.get(key) != value for key, value in expected.items()):
        raise ValueError("Native phylogeny configuration differs from planned command")
    proteomes = provenance["input_proteomes"]
    taxa = {item["filename"]: item["taxon"] for item in proteomes}
    if len(taxa) != len(proteomes) or len(set(taxa.values())) != len(taxa) or set(taxa) != set(fetch_fasta_files(str(fasta))):
        raise ValueError("Native/FASTA species inventory differs")
    tree_taxa = [t.name for t in Phylo.read(phylo / "species_tree.rooted.nwk", "newick").get_terminals()]
    if len(tree_taxa) != len(taxa) or set(tree_taxa) != set(taxa.values()):
        raise ValueError("Species tree taxa differ from inputs")
    owners, checked = {}, [*outputs.values(), record(phylo / "provenance_manifest.json"),
        record(phylo / "reconciliation_summary.json"), record(phylo / "species_tree.rooted.nwk")]
    for item in proteomes:
        identity = record(fasta / item["filename"])
        if identity["sha256"] != item["sha256"]:
            raise ValueError("Native FASTA checksum differs")
        checked.append(identity)
        for sequence in SeqIO.parse(identity["path"], "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty input sequence")
            owners[sequence.id] = item["taxon"]
    candidate = Path(argument(argv, "--candidate-clusters"))
    prediction = phylo.parent / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    candidates = {g: i for i, group in enumerate(read_partition(candidate, set(owners))) for g in group}
    groups = read_partition(prediction, set(owners))
    roots = read_predictions(phylo / PRODUCTS[0], "root_hogs")
    if set(membership(roots)) != set(owners) or {frozenset(g) for g in groups} != {frozenset(g) for g in roots.values()}:
        raise ValueError("Raw partition and RootHOG outputs differ")
    summary = json.loads((phylo / "reconciliation_summary.json").read_text())
    pairs = check_pairs(phylo / PRODUCTS[1], owners, candidates, summary["ortholog_pairs"])
    if pairs <= 0 or validate_confidence(phylo / PRODUCTS[2], phylo / PRODUCTS[1]) != pairs:
        raise ValueError("Invalid native pair/confidence counts")
    if summary["root_hogs"] != len(groups) or summary["candidate_families"] != len(set(candidates.values())):
        raise ValueError("Native family counts differ")
    if (row["coverage"].get("all_candidate_genes_preserved") is not True
            or len(owners) != row["coverage"]["candidate_genes"] or len(owners) != row["coverage"]["prediction_genes"]):
        raise ValueError("Claimed gene coverage differs")
    checked.extend([record(candidate), record(prediction)])
    return {"genes": len(owners), "species": len(taxa), "groups": len(groups), "pairs": pairs,
            "phylogeny": str(phylo), "prediction": str(prediction), "products": outputs,
            "checked": checked, "summary": summary}


def readback_scores(repo, dataset, row, native, folder, participant):
    checked = list(row.get("references", []))
    for item in checked:
        check(item)
    if dataset == "OrthoBench":
        result = score_orthobench(repo, Path(native["prediction"]))
        checked.extend(result["references"])
        details = {"precision": result["details"]["precision"] / 100,
                   "recall": result["details"]["recall"] / 100, "family_details": result["details"]}
        scores = result["scores"]
    elif dataset == "ThreeKingdoms":
        score_path = folder / "busco_score.txt"
        reported = parse_score(score_path)
        checked.append(record(score_path))
        script = repo / "three_kingdoms/score_against_busco.py"
        reference = repo / "three_kingdoms/busco/reference_orthogroups.txt"
        command = [sys.executable, str(script), "--predictions", native["prediction"], "--reference", str(reference)]
        text = subprocess.check_output(command, text=True)
        # Keep native result directories read-only; parse the independent run in memory.
        with TemporaryDirectory(prefix="ohmm-score-readback-") as temporary:
            repeated_path = Path(temporary) / "score.txt"
            repeated_path.write_text(text)
            repeated = parse_score(repeated_path)
        if repeated != reported:
            raise ValueError("Three Kingdoms score does not reproduce")
        checked.extend([record(script), record(reference)])
        details, scores = repeated, {"Three Kingdoms F1": repeated["f_score"]}
    else:
        conversion = row["conversion"]
        if conversion != {"written_pairs": native["pairs"], "total_pairs": native["pairs"],
                          "retained_pairs": native["pairs"], "removed_mapping_pairs": 0}:
            raise ValueError("QfO conversion count/mapping loss differs")
        for item in row["metric_files"]:
            check(item)
        checked.extend(row["metric_files"])
        reference = repo / "qfo_benchmark/benchmark-webservice/reference_data"
        assessment, paths = validate_directory(folder / "qfo_scoring", participant, reference)
        checked.extend(record(path) for path in paths)
        scores, details = {}, {}
        for challenge, result in assessment["endpoints"].items():
            metric = METRIC_NAMES["QfO"][challenge]
            scores[metric] = result["score"]
            axes = result["native_participant"]
            details[metric] = ({"precision": axes["metric_y"], "recall": axes["metric_x"]} if metric.endswith(" F")
                              else {"assessed_relations": axes["metric_x"], "similarity": axes["metric_y"]})
    if scores != row["scores"]:
        raise ValueError("Independent scores disagree with claimed scores")
    endpoint_vector(dataset, scores)
    return scores, details, checked


def audit_cell(repo, kind, dataset, row, execution):
    if row["status"] != "scored_development_result":
        raise ValueError("Cannot admit an unscored cell")
    verified = "inputs_sources_runtime_reverified" if kind == "full_search" else "inputs_reverified_unchanged"
    if row.get(verified) is not True:
        raise ValueError("Missing producer postflight validation")
    checked = [*execution["sources"]]
    if kind == "supplement":
        case = {"command": row["native"]["command"]}
        folder = Path(argument(case["command"], "--json")).parent
        checked.extend([*execution["inputs"], *execution["runtime"]["files"]])
        native = replay_native(row, case)
        if row["label"] == "control_direct":
            compare_golden(Path(native["phylogeny"]), execution["golden"])
            for key, item in (("candidate_partition", execution["baseline_candidate"]),
                              ("membership_constraints", execution["baseline_constraints"])):
                flag = "--candidate-clusters" if key == "candidate_partition" else "--membership-constraints"
                if record(argument(case["command"], flag))["sha256"] != item["sha256"]:
                    raise ValueError("Supplement candidate control differs")
        participant = None
    else:
        cases = [c for c in execution["cases"] if c["dataset"] == dataset and c["label"] == row["label"]]
        if len(cases) != 1:
            raise ValueError("Ambiguous experiment cell")
        case, folder = cases[0], Path(cases[0]["output"])
        if kind == "full_search":
            from benchmark_tools.run_full_search_hillclimb import checked_records
            checked.extend(checked_records(execution, case))
            result = validate_native(case, execution["cpu"])
            phylo = folder / "output/orthohmm_phylogeny"
            native = {**result, "phylogeny": str(phylo), "prediction": str(folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"),
                      "checked": list(result["products"].values())}
            participant = "ohmm_hillclimb_r1_" + case["scoring_id"]
        else:
            checked.extend([execution["plan"], *case["inputs"], *case["golden"].values(), *case["scoring_inputs"]])
            native = replay_native(row, case)
            if case["membership_rule"] == "high_confidence_pair":
                compare_golden(Path(native["phylogeny"]), case["golden"])
            participant = "ohmm_hillclimb_r1_" + str(case["index"])
    if native["genes"] != EXPECTED_GENES[dataset]:
        raise ValueError("Unexpected benchmark gene universe")
    checked.extend(native["checked"])
    for item in checked:
        check(item)
    scores, details, scoring_checked = readback_scores(repo, dataset, row, native, folder, participant)
    checked.extend(scoring_checked)
    for item in checked:
        check(item)
    return {"status": "independent_content_readback_verified", "scores": endpoint_vector(dataset, scores),
        "details": details, "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
        "products": native["products"], "checked": checked,
        "resources": {key: row[key] for key in ("resources", "native_metrics", "graph", "candidate_resources",
                      "total_incremental_wall_s", "total_full_pipeline_and_scoring_wall_s") if key in row},
        "scope": "Content and score readback; producer source/input pins rechecked. Not terminal scheduler admission or universal runtime closure."}


def combine_cells(methods, cells):
    for cell in cells:
        key = cell["method_key"]
        method = methods.setdefault(key, {"key": key, "label": key, "scores": {}, "cells": []})
        method["cells"].append(cell)
        if cell["status"] != "independent_content_readback_verified":
            continue
        existing = [c for c in method["cells"][:-1] if c["dataset"] == cell["dataset"]
                    and c["status"] == "independent_content_readback_verified"]
        if any(c["scores"] != cell["scores"] or {k: v["sha256"] for k, v in c["products"].items() if k in PRODUCTS}
               != {k: v["sha256"] for k, v in cell["products"].items() if k in PRODUCTS} for c in existing):
            raise ValueError("Conflicting original/supplemental metric cell")
        method["scores"].update(cell["scores"])
    return methods


def screen(vectors, control):
    baseline = vectors.get(control, {})
    comparisons = {}
    for label, scores in vectors.items():
        deltas = {key: scores[key] - baseline[key] for key in scores.keys() & baseline.keys()}
        complete = set(scores) == set(baseline) == set(ENDPOINTS)
        comparisons[label] = {"deltas_from_matched_control": deltas, "complete_eight_endpoint_contrast": complete,
            "passes_point_estimate_screen": label != control and complete
                and any(value > 0 for key, value in deltas.items() if key in PRIMARY_ENDPOINTS)
                and all(value >= -.005 for value in deltas.values())}
    return {"control": control, "comparisons": comparisons, "pareto": metric_frontier(vectors),
            "scope": "Development point-estimate screening only, not default promotion, statistical equivalence or independent confirmation."}


def build_report(repo, registry_path, membership_root, supplement_path, full_root, output,
                 deferred_fixed_root=None, deferred_inferred_root=None, deferred_qfo_root=None, deferred_tk_root=None):
    output.mkdir(parents=True, exist_ok=False)
    snapshots = output / "snapshots"
    snapshots.mkdir()
    sources = reporter_sources()
    registry, identity = snapshot(registry_path, snapshots, "retained_registry.json")
    methods, captures, cells = retained_methods(registry), [identity], []
    specs = [("membership", dataset, membership_root / dataset / "results.json") for dataset in DATASETS]
    specs += [("supplement", "ThreeKingdoms", supplement_path)]
    specs += [("full_search", dataset, full_root / dataset / "results.json") for dataset in DATASETS]
    for kind, dataset, path in specs:
        if not path.exists():
            for label in ROUND_ARMS["membership" if kind == "supplement" else kind]:
                if kind == "supplement" and label == "control_path":
                    continue
                cells.append({"method_key": ("membership" if kind == "supplement" else kind) + ":" + label,
                    "dataset": dataset, "status": "report_missing_execution_state_not_inferred", "scores": {}})
            continue
        report, capture = snapshot(path, snapshots, f"{kind}_{dataset}.json")
        captures.append(capture)
        execution = None
        try:
            check(report["execution_plan"])
            execution, plan_capture = snapshot(report["execution_plan"]["path"], snapshots, f"{kind}_{dataset}_plan.json")
            captures.append(plan_capture)
            if execution["source_revision"] != report["source_revision"]:
                raise ValueError("Producer/execution revision differs")
        except Exception as error:
            execution_error = str(error)
        observed = set()
        for row in report["arms"]:
            label = row["label"]
            if label in observed or label not in ROUND_ARMS["membership" if kind == "supplement" else kind]:
                raise ValueError("Duplicate or unknown arm label")
            observed.add(label)
            cell = {"method_key": ("membership" if kind == "supplement" else kind) + ":" + label,
                    "dataset": dataset, "producer_status": row["status"], "status": row["status"], "scores": {},
                    "source_snapshot": capture["snapshot"], "claimed_scores": row.get("scores", {}),
                    "producer_job_id": report.get("job_id"), "producer_scope": report.get("scope")}
            if row["status"] == "scored_development_result":
                write_json(output / "progress.json", {"status": "content_readback_running",
                    "active_cell": {"kind": kind, "dataset": dataset, "label": label}, "processed_cells": len(cells)})
                try:
                    if execution is None:
                        raise ValueError(execution_error)
                    cell.update(audit_cell(repo, kind, dataset, row, execution))
                except Exception as error:
                    cell.update(status="readback_failed", error_type=type(error).__name__, error=str(error))
            cells.append(cell)
        for label in ROUND_ARMS["membership" if kind == "supplement" else kind]:
            if label not in observed and not (kind == "supplement" and label == "control_path"):
                cells.append({"method_key": ("membership" if kind == "supplement" else kind) + ":" + label,
                    "dataset": dataset, "status": "not_yet_reported_execution_state_not_inferred", "scores": {}})
    deferred_specs = []
    if deferred_fixed_root is not None:
        deferred_specs.append(("deferred_fixed", "QfO", deferred_fixed_root))
    if deferred_inferred_root is not None:
        deferred_specs.extend(("deferred_inferred", dataset,
                               deferred_tk_root if dataset == "ThreeKingdoms" and deferred_tk_root else deferred_inferred_root)
                              for dataset in ("OrthoBench", "ThreeKingdoms"))
    if deferred_qfo_root is not None:
        deferred_specs.append(("deferred_inferred", "QfO", deferred_qfo_root))
    for kind, dataset, root in deferred_specs:
        from benchmark_tools.report_deferred_cells import collect_deferred_cells
        write_json(output / "progress.json", {"status": "content_readback_running",
            "active_cell": {"kind": kind, "dataset": dataset}, "processed_cells": len(cells)})
        admitted, extra_captures = collect_deferred_cells(repo, kind, dataset, root, snapshots)
        cells.extend(admitted)
        captures.extend(extra_captures)
    combine_cells(methods, cells)
    for item in sources:
        check(item)
    screenings = {}
    rounds = dict(ROUND_ARMS)
    rounds.update({kind: ("control", "deferred500") for kind, _, _ in deferred_specs})
    for kind, labels in rounds.items():
        vectors = {label: methods.get(kind + ":" + label, {"scores": {}})["scores"] for label in labels}
        screenings[kind] = screen(vectors, labels[0])
    report = {"schema_version": 1, "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "reporter_sources": sources,
        "status": "partial_or_complete_vectors_described_not_goal_completion", "captures": captures,
        "methods": list(methods.values()), "round_screening": screenings, "defaults_promoted": False,
        "readback_failures": [c for c in cells if c["status"] == "readback_failed"],
        "limitations": [*registry["limitations"], "All evaluated datasets are development-exposed.",
            "TreeFam is supporting evidence. Functional similarity is not orthology correctness; Three Kingdoms is BUSCO-reference co-membership.",
            "FAS is unseeded: small point differences and Pareto dominance involving it are descriptive, not demonstrated superiority.",
            "Resource observations retain their original scopes; cached replay and fresh full-pipeline times must not be ranked together.",
            "Fixed-tree copy-deferral diagnostics are separate from candidate-specific inferred-tree vectors and cannot fill their missing endpoints.",
            "Incomplete cells remain missing, not zero. Content snapshots do not imply all scheduler tasks have completed."]}
    write_json(output / "report.json", report)
    write_json(output / "progress.json", {"status": "content_snapshot_readback_finished",
        "processed_cells": len(cells), "readback_failures": len(report["readback_failures"]), "goal_completed": False})
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["method", *ENDPOINTS, "complete_eight_endpoints"])
        for method in methods.values():
            writer.writerow([method["label"], *(method["scores"].get(key, "NA") for key in ENDPOINTS),
                             set(method["scores"]) == set(ENDPOINTS)])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "registry", "membership-root", "supplement", "full-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("deferred-fixed-root", "deferred-inferred-root", "deferred-qfo-root", "deferred-tk-root"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    report = build_report(args.repo.resolve(), args.registry.resolve(), args.membership_root.resolve(),
                          args.supplement.resolve(), args.full_root.resolve(), args.output.resolve(),
                          args.deferred_fixed_root.resolve() if args.deferred_fixed_root else None,
                          args.deferred_inferred_root.resolve() if args.deferred_inferred_root else None,
                          args.deferred_qfo_root.resolve() if args.deferred_qfo_root else None,
                          args.deferred_tk_root.resolve() if args.deferred_tk_root else None)
    print(json.dumps({"methods": len(report["methods"]), "readback_failures": len(report["readback_failures"])}))
    return 1 if report["readback_failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
