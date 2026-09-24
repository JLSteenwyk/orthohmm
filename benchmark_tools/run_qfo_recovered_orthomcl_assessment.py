"""Score admitted recovered OrthoMCL group-clique pairs on frozen QfO endpoints."""

import argparse
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_orthomcl import completed, frozen, unique_records
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_recovered_orthomcl_pairs import validate_admission, SEMANTICS, ADMITTER_SHA
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen

CONVERTER_SHA = "ceaf8898c373aa9c0ee75399682c3bd8a509f09e13c012c1623238f438e7a5fc"


def validate_stage(stage, scheduler):
    if (stage["status"] != "recovered_orthomcl_pairs_prepared_unscored"
            or stage["accuracy_evaluated"] is not False or stage["publication_ready"] is not False
            or stage["method"] != "orthomcl" or stage["participant"] != "qfo_corrected_orthomcl_recovered"
            or stage["semantics"] != SEMANTICS or stage["job_id"] != scheduler["JobIDRaw"]):
        raise ValueError("Wrong recovered conversion identity/status")
    expected = dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="64G")
    if any(scheduler.get(k) != v for k, v in expected.items()):
        raise ValueError("Require completed recovered pair conversion")
    for key in ("total_pairs", "retained_pairs", "removed_mapping_pairs", "native_duplicate_relations"):
        if type(stage[key]) is not int or stage[key] < 0:
            raise ValueError("Invalid recovered pair accounting")
    if (stage["total_pairs"] <= 0 or stage["retained_pairs"] != stage["total_pairs"]
            or stage["removed_mapping_pairs"] != 0 or stage["native_duplicate_relations"] != 0
            or any(stage["pairs"][k] != stage["filtered_pairs"][k] for k in ("bytes", "sha256"))):
        raise ValueError("Changed/lost recovered reference pairs")
    content = stage["content"]
    if (set(content) != {"total_pairs", "final_groups", "grouped_proteins", "ungrouped_input_proteins"}
            or any(type(v) is not int or v < 0 for v in content.values())
            or content["total_pairs"] != stage["total_pairs"] or content["final_groups"] <= 0
            or content["grouped_proteins"] + content["ungrouped_input_proteins"] != 984137):
        raise ValueError("Invalid recovered clique/coverage accounting")


def prepare(root, digest, job, executor, commit, require_fresh=True, helpers=None):
    helpers = Path(__file__).resolve().parent if helpers is None else helpers
    scheduler, accounting = completed(job, 2, "64G")
    if not executor.resolve().is_relative_to(root / "benchmarks/work"):
        raise ValueError("Require retained conversion executor")
    frozen(executor, commit)
    directory = root / "benchmarks/results/qfo_blast_recovery_pairs_v1"
    path = directory / "results.json"
    stage = read_frozen(path, digest)
    validate_stage(stage, scheduler)
    source = record(executor / "benchmark_tools/prepare_recovered_orthomcl_pairs.py")
    if source["sha256"] != CONVERTER_SHA or stage["source"] != source or source not in stage["checked_records"]:
        raise ValueError("Unbound recovered converter source")
    for key, expected in (("pairs", directory / "pairs.tsv"), ("filtered_pairs", directory / "pairs.qfo.tsv"),
        ("group_audit", directory / "groups.json"),
        ("admission", root / "benchmarks/results/qfo_blast_recovery_native_admission_v1/report.json")):
        if stage[key]["path"] != str(expected):
            raise ValueError("Unexpected recovered pair artifact path")
        check(stage[key])
    if stage["admission"] not in stage["checked_records"]:
        raise ValueError("Unbound recovered native admission")
    prior = stage["admission"]
    admission = read_frozen(Path(prior["path"]), prior["sha256"])
    prior_scheduler, _ = completed(int(admission["admission_job_id"]), 2, "64G")
    content = validate_admission(admission, admission["admission_job_id"])
    if (prior_scheduler != stage["admission_scheduler"] or admission["source"]["sha256"] != ADMITTER_SHA
            or admission["query_coverage"] != stage["query_coverage"]):
        raise ValueError("Changed recovered native admission or failure coverage")
    groups = read_frozen(Path(stage["group_audit"]["path"]), stage["group_audit"]["sha256"])
    if (groups["status"] != "native_final_groups_match_mcl_partition" or groups["content"] != content
            or groups["accuracy_admitted"] is not False or groups["publication_ready"] is not False):
        raise ValueError("Changed recovered final-group audit")
    expected = {k: content[k] for k in ("final_groups", "grouped_proteins", "ungrouped_input_proteins")}
    expected["total_pairs"] = content["cross_species_clique_pairs"]
    if expected != stage["content"]:
        raise ValueError("Recovered conversion differs from native audit")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    manifest = read_frozen(env_path, ENV_SHA)
    if [r for r in manifest["reference_files"] if Path(r["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
        raise ValueError("Conversion/scoring reference mapping differs")
    records = unique_records([record(path), record(env_path), *environment_records(manifest),
        source, *stage["checked_records"], stage["pairs"], stage["filtered_pairs"], stage["group_audit"],
        *groups["checked_records"], admission["source"], *admission["checked_records"], *admission["outputs"],
        *[record(p) for p in sorted(helpers.glob("*.py"))]])
    for item in records:
        check(item)
    output = root / "benchmarks/results/qfo_blast_recovery_assessment_v1"
    work = root / "qfo_benchmark/w/qc_mcr"
    results = root / "qfo_benchmark/scoring/corrected_orthomcl_recovered"
    for destination in (output, work, results):
        if require_fresh and (destination.exists() or destination.is_symlink()):
            raise FileExistsError(destination)
    return dict(status="prepared_unrun", method="orthomcl", stage=stage,
        source=record(helpers / "run_qfo_recovered_orthomcl_assessment.py"),
        pairs_manifest=record(path), environment_manifest=record(env_path), conversion_scheduler=scheduler,
        conversion_accounting=accounting, command=command_for(root, stage, manifest, work, results),
        cwd=str(output), work=str(work), results=str(results), verified_records=records,
        environment_overrides=manifest["environment_overrides"], accuracy_admitted=False, publication_ready=False)


def run(root, digest, job, executor, commit, check_only=False):
    report = prepare(root, digest, job, executor, commit)
    if check_only:
        return report
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "8"
            or os.environ.get("SLURM_MEM_PER_NODE") != "65536" or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled eight-CPU/64-GiB assessment on bizon")
    output = Path(report["cwd"])
    output.mkdir(parents=True, exist_ok=False)
    report.update(status="running", job_id=os.environ["SLURM_JOB_ID"])
    save_status(output / "preflight.json", report)
    env = {**os.environ, **report["environment_overrides"],
           "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    try:
        with (output / "scoring.log").open("x") as log:
            done = subprocess.run(report["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = done.returncode
        for item in [report["source"], *report["verified_records"]]:
            check(item)
        report["outputs"] = [record(p) for p in sorted(Path(report["results"]).rglob("*")) if p.is_file()]
        if done.returncode:
            raise RuntimeError(f"Recovered native scoring failed: {done.returncode}")
        report["status"] = "process_succeeded_pending_independent_admission"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        if (output / "scoring.log").exists():
            report["log"] = record(output / "scoring.log")
        save_status(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "executor"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--conversion-job", type=int, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    print(run(args.root.resolve(), args.pairs_sha256, args.conversion_job,
              args.executor.resolve(), args.commit, args.check_only)["status"])
