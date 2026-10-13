"""Prepare and execute the bounded candidate/membership development factorial."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.metric_frontier import ENDPOINTS, PRIMARY_ENDPOINTS, metric_frontier

DATASETS = ("OrthoBench", "QfO", "ThreeKingdoms")
PRODUCTS = ("orthohmm_root_hogs.tsv", "orthohmm_pairwise_orthologs.tsv",
            "orthohmm_pairwise_orthologs_confidence.tsv")


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(temporary, path)


def argument(argv, flag):
    if argv.count(flag) != 1:
        raise ValueError(f"Require exactly one {flag}")
    return argv[argv.index(flag) + 1]


def derive_command(original, source, output, policy, cpu):
    argv = list(original)
    if argument(argv, "--species-tree-mode") != "infer":
        raise ValueError("Require independently inferred baseline species tree")
    argv[:2] = [sys.executable, str(source / "benchmark_tools/replay_phylogeny.py")]
    for flag, value in (("--output-directory", output / "output"),
                        ("--json", output / "metrics.json"), ("--cpu", cpu)):
        argv[argv.index(flag) + 1] = str(value)
    if "--checkpoint-source" in argv:
        index = argv.index("--checkpoint-source")
        del argv[index:index + 2]
    # Each candidate arm reuses its own completed native trees, not another partition.
    argv.extend(["--checkpoint-source", argument(original, "--output-directory"),
                 "--satellite-membership-rule", policy])
    return argv


def original_commands(repo, dataset):
    results = repo / "benchmark_tools/results"
    if dataset == "OrthoBench":
        path = results / "orthobench_factorial_prepared_20260916.json"
        prepared = json.loads(path.read_text())
        control = next(c["argv"] for c in prepared["cells"] if c["label"] == "p1_c1_r1")
        low_path = repo / "benchmarks/results/ob_cpm_phylogeny_v1/cpm_low/preflight.json"
        broader = json.loads(low_path.read_text())["cell"]["argv"]
        return control, broader, [record(path), record(low_path)]
    if dataset == "QfO":
        path = repo / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
        prepared = json.loads(path.read_text())
        control = next(c["argv"] for c in prepared["cells"] if c["label"] == "p1_c1_r1")
        low_path = repo / "benchmarks/results/qfo_cpm_phylogeny_v1/cpm_low/preflight.json"
        broader = json.loads(low_path.read_text())["executed_argv"]
        return control, broader, [record(path), record(low_path)]
    if dataset != "ThreeKingdoms":
        raise ValueError("Unknown dataset")
    base = repo / "three_kingdoms/results/parity_20260907/orthohmm_phylogeny_satellite_v2"
    metrics_path = base / "metrics.json"
    metadata = json.loads(metrics_path.read_text())["metadata"]
    output = Path(metadata["output_directory"])
    control = [sys.executable, str(repo / "benchmark_tools/replay_phylogeny.py"),
        "--fasta-directory", metadata["fasta_directory"],
        "--candidate-clusters", str(output / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"),
        "--membership-constraints", str(output / "orthohmm_working_res/phylogeny_candidate_merges.json"),
        "--output-directory", str(output), "--json", str(metrics_path), "--cpu", "32",
        "--species-tree-mode", "infer", "--aligner", metadata["aligner"],
        "--tree-builder", metadata["tree_builder"], "--root-rule", "species_overlap",
        "--pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance"]
    return control, None, [record(metrics_path)]


def scoring_records(repo, dataset):
    if dataset == "OrthoBench":
        benchmark = repo.parent.parent / "Open_Orthobench/BENCHMARKS"
        return [record(benchmark / "benchmark.py"),
                *[record(p) for p in sorted((benchmark / "RefOGs").rglob("*.txt"))]]
    if dataset == "ThreeKingdoms":
        return [record(repo / "three_kingdoms/score_against_busco.py"),
                record(repo / "three_kingdoms/busco/reference_orthogroups.txt")]
    from benchmark_tools.run_qfo_recovered_assessment import environment_records
    path = repo / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    return [record(path), *environment_records(json.loads(path.read_text()))]


def prepare(repo, source, output, cpu):
    if cpu < 1:
        raise ValueError("CPU budget must be positive")
    if output.exists():
        raise FileExistsError(output)
    plan_path = source / "benchmark_tools/hillclimb_round1_20261012.json"
    plan = json.loads(plan_path.read_text())
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools/run_membership_hillclimb.py",
                    "benchmark_tools/replay_phylogeny.py"], check=True)
    cases = []
    for dataset in DATASETS:
        control, broader, evidence = original_commands(repo, dataset)
        scoring_inputs = scoring_records(repo, dataset)
        for index, arm in enumerate(plan["arms"]):
            old = control if arm["cpm_resolution"] == .1 else broader
            row = {"dataset": dataset, "index": index, **arm}
            if old is None:
                row.update(status="unavailable", reason="No completed CPM 0.08 Three Kingdoms candidate/tree checkpoint. A separate preparation is required.")
                cases.append(row)
                continue
            folder = output / dataset / arm["label"]
            argv = derive_command(old, source, folder, arm["membership_rule"], cpu)
            tools = []
            for flag in ("--aligner", "--tree-builder"):
                executable = shutil.which(argument(argv, flag))
                if executable is None:
                    raise FileNotFoundError(argument(argv, flag))
                # Keep original tool spelling: it participates in cache hashes.
                tools.append(record(executable))
            old_output = Path(argument(old, "--output-directory"))
            golden = {name: record(old_output / "orthohmm_phylogeny" / name) for name in PRODUCTS}
            golden["species_tree.rooted.nwk"] = record(old_output / "orthohmm_phylogeny/species_tree.rooted.nwk")
            inputs = [*evidence, *tools, record(Path(argument(old, "--candidate-clusters"))),
                      record(Path(argument(old, "--membership-constraints"))),
                      record(old_output / "orthohmm_phylogeny/provenance_manifest.json"),
                      *[record(p) for p in sorted(Path(argument(old, "--fasta-directory")).iterdir()) if p.is_file()]]
            row.update(status="prepared", command=argv, output=str(folder), inputs=inputs,
                       golden=golden, historical_command=old, scoring_inputs=scoring_inputs)
            cases.append(row)
    source_files = [*[source / "benchmark_tools" / name for name in (
                        "run_membership_hillclimb.py", "replay_phylogeny.py",
                        "prepare_ob_candidate_neighborhood.py", "orthobench_stage_diagnostics.py",
                        "score_orthobench_partition.py", "summarize_three_kingdoms_parity.py",
                        "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py", "qfo_summarize_scores.py",
                        "run_qfo_recovered_assessment.py", "simulation_method_outputs.py")],
                    source / "qfo_benchmark/og_to_pairwise.py",
                    *sorted((source / "orthohmm").glob("*.py"))]
    report = {"schema_version": 1, "status": "prepared_unrun", "plan": record(plan_path),
              "source_root": str(source), "source_revision": revision,
              "sources": [record(p) for p in source_files], "repo": str(repo),
              "output": str(output), "cpu": cpu, "cases": cases,
              "tool_path": os.environ["PATH"],
              "scope": "Development-exposed cached-HMM/candidate and native reconciliation replay; not fresh initial-search recovery or isolated full-pipeline timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", report)
    return report


def compare_golden(phylogeny, golden):
    matches = {name: record(phylogeny / name)["sha256"] == item["sha256"]
               for name, item in golden.items()}
    if not all(matches.values()):
        raise ValueError("Direct control failed exact retained-output check: " + str(matches))
    return matches


def validate_partition(candidate, prediction):
    def genes(path):
        seen = set()
        for line in path.read_text().splitlines():
            for gene in line.split():
                if gene in seen:
                    raise ValueError(f"Duplicate partition gene: {gene}")
                seen.add(gene)
        return seen
    original, inferred = genes(candidate), genes(prediction)
    if original != inferred:
        raise ValueError("Candidate/prediction gene universes differ")
    return {"candidate_genes": len(original), "prediction_genes": len(inferred),
            "all_candidate_genes_preserved": True}


def score_orthobench(repo, prediction):
    from benchmark_tools.orthobench_stage_diagnostics import read_clusters, run_official_benchmark
    from benchmark_tools.score_orthobench_partition import read_named_groups, score_partition
    benchmark = repo.parent.parent / "Open_Orthobench/BENCHMARKS"
    refogs = benchmark / "RefOGs"
    names = sorted(p.name for p in refogs.glob("RefOG*.txt"))
    if not names:
        raise ValueError("No OrthoBench reference groups")
    references = read_named_groups(refogs, names)
    uncertain = {name: set((refogs / "low_certainty_assignments" / name).read_text().split())
                 if (refogs / "low_certainty_assignments" / name).is_file() else set() for name in names}
    score = score_partition(read_clusters(prediction), references, uncertain)
    official = run_official_benchmark(benchmark / "benchmark.py", prediction)
    # The official script prints percentages to one decimal place.
    if any(abs(official[key] - score[key]) > .05000001 for key in ("f_score", "precision", "recall")):
        raise ValueError("Independent and official OrthoBench scores disagree")
    return {"scores": {"OrthoBench F1": score["f_score"] / 100},
            "details": score, "official": official,
            "references": [record(benchmark / "benchmark.py"), *[record(refogs / name) for name in names],
                           *[record(p) for p in sorted((refogs / "low_certainty_assignments").glob("*.txt"))]]}


def score_three_kingdoms(repo, prediction, folder):
    from benchmark_tools.summarize_three_kingdoms_parity import parse_score
    reference = repo / "three_kingdoms/busco/reference_orthogroups.txt"
    script = repo / "three_kingdoms/score_against_busco.py"
    path = folder / "busco_score.txt"
    command = [sys.executable, str(script), "--predictions", str(prediction), "--reference", str(reference)]
    with path.open("x") as stream:
        subprocess.run(command, check=True, stdout=stream, text=True)
    score = parse_score(path)
    return {"scores": {"Three Kingdoms F1": score["f_score"]}, "details": score,
            "references": [record(reference), record(script)], "command": command}


def score_qfo(repo, phylogeny, folder, index):
    from Bio import SeqIO
    from benchmark_tools.prepare_qfo_factorial_pairs import write_native_pairs
    from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
    from benchmark_tools.qfo_summarize_scores import summarize, METRICS, metric_path
    from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records, require_darwin_work_path
    work = repo / "qfo_benchmark/w" / f"hc1_{index}"
    results = folder / "qfo_scoring"
    require_darwin_work_path(work)
    if work.exists() or results.exists():
        raise FileExistsError("QfO scoring namespace already exists")
    environment_path = repo / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = json.loads(environment_path.read_text())
    identities = [record(environment_path), *environment_records(environment)]
    for item in identities:
        check(item)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Ambiguous QfO mapping")
    manifest = json.loads((phylogeny / "provenance_manifest.json").read_text())
    fasta = repo / "benchmarks/work/qfo_corrected_inputs_20260918"
    owners = {}
    for item in manifest["input_proteomes"]:
        for protein in SeqIO.parse(fasta / item["filename"], "fasta"):
            if protein.id in owners:
                raise ValueError("Duplicate QfO input gene")
            owners[protein.id] = item["taxon"]
    count = json.loads((phylogeny / "reconciliation_summary.json").read_text())["ortholog_pairs"]
    pairs, filtered = folder / "pairs.tsv", folder / "pairs.qfo.tsv"
    written = write_native_pairs(phylogeny / "orthohmm_pairwise_orthologs.tsv", pairs, owners, count)
    total, retained = filter_pairs(pairs, filtered, load_mapping(Path(mappings[0]["path"])))
    if not 0 < written == total == retained:
        raise ValueError("Corrected QfO pair conversion has mapping loss")
    stage = {"filtered_pairs": record(filtered), "participant": f"ohmm_hillclimb_r1_{index}"}
    command = command_for(repo, stage, environment, work, results)
    env = {**os.environ, **environment["environment_overrides"],
           "NXF_SINGULARITY_CACHEDIR": str(repo / "qfo_benchmark/scoring/container_cache")}
    write_json(folder / "scoring_preflight.json", {"command": command, "inputs": identities,
               "stage": stage, "total_pairs": total, "retained_pairs": retained})
    with (folder / "qfo_scoring.log").open("x") as log:
        subprocess.run(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    for item in identities:
        check(item)
    values = summarize(results)
    if any(not math.isfinite(value) or not 0 <= value <= 1 for value in values["scores"].values()):
        raise ValueError("Invalid QfO endpoint")
    return {"scores": values["scores"], "metric_files": [record(metric_path(results, path)) for path, _ in METRICS.values()],
            "conversion": {"written_pairs": written, "total_pairs": total, "retained_pairs": retained,
                           "removed_mapping_pairs": total - retained}, "command": command,
            "scope": "Native descriptive endpoints; FAS sampling is unseeded. Not an independent scoring admission."}


def summarize_round(output, summary_name="summary.json"):
    vectors = {}
    states = []
    for dataset in DATASETS:
        path = output / dataset / "results.json"
        if not path.exists():
            states.append({"dataset": dataset, "status": "pending"})
            continue
        report = json.loads(path.read_text())
        states.append({"dataset": dataset, "status": report["status"], "receipt": record(path)})
        for row in report["arms"]:
            if row.get("status") == "scored_development_result":
                values = vectors.setdefault(row["label"], {})
                for metric, value in row["scores"].items():
                    values[dataset + ": " + metric] = value
    baseline = vectors.get("control_direct", {})
    labels = ("control_direct", "control_path", "broader_direct", "broader_path")
    comparisons = []
    for label in labels:
        values = vectors.get(label, {})
        deltas = {metric: values[metric] - baseline[metric] for metric in values.keys() & baseline.keys()}
        complete = set(values) == set(baseline) == set(ENDPOINTS)
        primary = [v for key, v in deltas.items() if key in PRIMARY_ENDPOINTS]
        eligible = label != "control_direct" and complete and any(v > 0 for v in primary) and all(v >= -.005 for v in deltas.values())
        comparisons.append({"label": label, "scores": values, "deltas_from_control": deltas,
                            "all_eight_endpoints_available": complete,
                            "passes_point_estimate_development_screen": eligible})
    summary = {"status": "complete_metric_vectors" if all(r["all_eight_endpoints_available"] for r in comparisons) else "partial_metric_vectors",
               "datasets": states, "comparisons": comparisons, "defaults_promoted": False,
               "pareto": metric_frontier({label: vectors.get(label, {}) for label in labels}),
               "limitations": ["Development-exposed point estimates; screening is not statistical noninferiority or independent validation.",
                               "FAS sampling is unseeded; small differences may reflect sampling.",
                               "Shared-host incremental replay resource use is not full-pipeline efficiency evidence."]}
    write_json(output / summary_name, summary)
    return summary


def run(plan_path, dataset):
    plan = json.loads(plan_path.read_text())
    repo, source = Path(plan["repo"]), Path(plan["source_root"])
    if Path(__file__).resolve() != source / "benchmark_tools/run_membership_hillclimb.py":
        raise ValueError("Execute the frozen runner source")
    for item in [plan["plan"], *plan["sources"]]:
        check(item)
    folder = Path(plan["output"]) / dataset
    folder.mkdir(exist_ok=False)
    report = {"status": "running", "execution_plan": record(plan_path), "dataset": dataset,
              "job_id": os.environ.get("SLURM_JOB_ID"), "arms": [], "all_metric_screening_complete": False,
              "defaults_promoted": False, "source_revision": plan["source_revision"], "scope": plan["scope"]}
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    direct_valid = {}
    for case in [c for c in plan["cases"] if c["dataset"] == dataset]:
        row = {"label": case["label"], "membership_rule": case["membership_rule"],
               "cpm_resolution": case["cpm_resolution"], "scores": {}, "status": "running"}
        report["arms"].append(row)
        write_json(folder / "results.json", report)
        if case["status"] == "unavailable":
            row.update(status="unavailable", reason=case["reason"])
            write_json(folder / "results.json", report)
            continue
        arm_folder = Path(case["output"])
        started = time.monotonic()
        try:
            for item in [*case["inputs"], *case["golden"].values(), *case["scoring_inputs"]]:
                check(item)
            if case["membership_rule"] == "high_confidence_path" and not direct_valid.get(case["cpm_resolution"]):
                raise ValueError("Matching direct control failed; do not interpret path variant")
            arm_folder.mkdir(exist_ok=False)
            write_json(arm_folder / "preflight.json", case)
            with (arm_folder / "inference.log").open("x") as log:
                process = subprocess.run(case["command"], cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT)
            row["inference_exit_code"] = process.returncode
            if process.returncode:
                raise RuntimeError(f"Native replay failed: {process.returncode}")
            phylogeny = arm_folder / "output/orthohmm_phylogeny"
            prediction = arm_folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
            row["coverage"] = validate_partition(Path(argument(case["command"], "--candidate-clusters")), prediction)
            if case["membership_rule"] == "high_confidence_pair":
                row["golden_matches"] = compare_golden(phylogeny, case["golden"])
                direct_valid[case["cpm_resolution"]] = True
            else:
                if record(phylogeny / "species_tree.rooted.nwk")["sha256"] != case["golden"]["species_tree.rooted.nwk"]["sha256"]:
                    raise ValueError("Membership-only change altered species tree")
            row["native_summary"] = json.loads((phylogeny / "reconciliation_summary.json").read_text())
            row["membership_audit"] = json.loads((phylogeny / "provenance_manifest.json").read_text())["membership_reconciliation"]
            row["resources"] = json.loads((arm_folder / "metrics.json").read_text())
            row["outputs"] = [record(phylogeny / name) for name in PRODUCTS]
            row["status"] = "native_outputs_verified_scoring_pending"
            write_json(folder / "results.json", report)
            scoring = (score_orthobench(repo, prediction) if dataset == "OrthoBench"
                       else score_three_kingdoms(repo, prediction, arm_folder) if dataset == "ThreeKingdoms"
                       else score_qfo(repo, phylogeny, arm_folder, case["index"]))
            row.update(status="scored_development_result", **scoring)
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        finally:
            row["total_incremental_wall_s"] = time.monotonic() - started
            try:
                for item in [*case["inputs"], *case["golden"].values(), *case["scoring_inputs"], *plan["sources"]]:
                    check(item)
                row["inputs_reverified_unchanged"] = True
            except Exception as error:
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
            write_json(folder / "results.json", report)
    report["status"] = "finished_all_arms_scored" if all(r["status"] == "scored_development_result" for r in report["arms"]) else "finished_with_explicit_gaps"
    write_json(folder / "results.json", report)
    summarize_round(Path(plan["output"]))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    preparation = sub.add_parser("prepare")
    preparation.add_argument("--repo", type=Path, required=True)
    preparation.add_argument("--source-root", type=Path, required=True)
    preparation.add_argument("--output", type=Path, required=True)
    preparation.add_argument("--cpu", type=int, default=32)
    execution = sub.add_parser("run")
    execution.add_argument("--plan", type=Path, required=True)
    execution.add_argument("--dataset", choices=DATASETS, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(), args.cpu)
              if args.mode == "prepare" else run(args.plan.resolve(), args.dataset))
    print(json.dumps({"status": result["status"], "dataset": result.get("dataset")}, sort_keys=True))


if __name__ == "__main__":
    main()
