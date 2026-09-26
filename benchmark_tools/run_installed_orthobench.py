"""One fresh full OrthoBench run of the setup-overlay installed package."""

import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

PINS = {
    "publication_frozen_overlay_install_20260926.json": "467036ee88d5e6f15d322aea4c94397591ef7ba521ac9f8ddde4b7cb783b97d2",
    "publication_mafft_build_20260926.json": "4c4e92a29c1dc4b27e4b47f4010e9ea9648c17beb02b5a00995a59535394ac76",
    "publication_fasttree_acquisition_20260926.json": "c553af49b95153434c6c09e4bf4dd6bf215d472a593c6f13103394785570ec68",
    "orthobench_factorial_prepared_20260916.json": "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382",
    "orthobench_factorial_results_20260916.json": "6a0d588b5cb47c60fc6bc8bae8aa0c83e5f2aadb11de970919d8c6527c387141",
}


def fasta_ids(paths):
    genes = set()
    for path in paths:
        with path.open() as stream:
            for line in stream:
                if line.startswith(">"):
                    fields = line[1:].split()
                    if not fields or fields[0] in genes:
                        raise ValueError("Empty or duplicate FASTA identifier")
                    genes.add(fields[0])
    return genes


def command(python, fasta, destination, mafft, fasttree):
    return [str(python), "-I", "-m", "orthohmm", str(fasta), "-o", str(destination),
        "-c", "32", "--threads_per_worker", "1", "--search_mode", "builtin", "--clustering", "leiden",
        "--accuracy_profile", "high_sensitivity", "--phylogeny", "reconcile", "--species_tree_mode", "infer",
        "--species_tree_rooting", "min_variance", "--phylogeny_candidates", "satellite_v2",
        "--phylogeny_root_rule", "species_overlap", "--phylogeny_pair_rule", "positive_paralogy",
        "--aligner", str(mafft), "--tree_builder", str(fasttree)]


def prepare(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    reports, checked = {}, [record(__file__), record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"))]
    for name, expected in PINS.items():
        path = repo / "benchmark_tools/results" / name
        item = record(path)
        if item["sha256"] != expected:
            raise ValueError("Changed publication evidence: " + name)
        checked.append(item)
        reports[name] = json.loads(path.read_text())
    installed = reports["publication_frozen_overlay_install_20260926.json"]
    mafft = reports["publication_mafft_build_20260926.json"]
    fasttree = reports["publication_fasttree_acquisition_20260926.json"]["installed_binary"]
    inputs = reports["orthobench_factorial_prepared_20260916.json"]["fasta_inputs"]
    baseline = reports["orthobench_factorial_results_20260916.json"]
    checked.extend([*installed["checked_records"], mafft["launcher"], *mafft["built_helpers"],
                    fasttree, *inputs, baseline["predictions"]["p1_c1_r1"], *baseline["references"]])
    for item in checked:
        check(item)
    if len(inputs) != 12 or len(fasta_ids([Path(x["path"]) for x in inputs])) != 251378:
        raise ValueError("Unexpected full OrthoBench universe")
    output.mkdir(parents=True)
    fasta = output / "input"
    fasta.mkdir()
    for item in inputs:
        source = Path(item["path"])
        destination = fasta / source.name
        shutil.copy2(source, destination)
        copied = record(destination)
        if (copied["bytes"], copied["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("FASTA copy differs")
        checked.append(copied)
    python = repo / "benchmarks/work/publication_frozen_overlay_20260926/venv_clean/bin/python"
    probe = json.loads(subprocess.check_output([str(python), "-I", "-c",
        "import json,sys,orthohmm;print(json.dumps({'module':orthohmm.__file__,'prefix':sys.prefix}))"],
        text=True, cwd=output))
    if not Path(probe["module"]).is_relative_to(python.parent.parent) or record(probe["module"]) not in checked:
        raise ValueError("Wrong installed import")
    plan = dict(status="prepared_not_run", checked_records=checked, imported_package=probe,
        command=command(python, fasta, output / "inference", mafft["launcher"]["path"], fasttree["path"]),
        cwd=str(output), executor=record(__file__),
        environment_overrides=dict(PATH="/usr/bin:/bin", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1", PYTHONHASHSEED="0", MAFFT_BINARIES=str(Path(mafft["launcher"]["path"]).parents[1] / "libexec/mafft")),
        environment_unset=["PYTHONPATH", "PYTHONHOME"],
        baseline_partition=baseline["predictions"]["p1_c1_r1"],
        baseline_scores=baseline["point_estimates_percent"]["p1_c1_r1"],
        expected_species=12, expected_genes=251378, attempts=1,
        resource_request=dict(cpus=32, memory_gib=128, time_hours=24),
        timing_scope="Shared-host fresh installation reproduction; not controlled comparative timing.",
        scientific_scores_admitted=False, publication_ready=False)
    with (output / "plan.json").open("x") as stream:
        json.dump(plan, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return plan


def run(plan_path):
    plan = json.loads(plan_path.read_text())
    if not os.environ.get("SLURM_JOB_ID") or int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) != 32:
        raise ValueError("Run only in the specified 32-CPU Slurm allocation")
    if record(__file__)["sha256"] != plan["executor"]["sha256"]:
        raise ValueError("Changed executor")
    for item in plan["checked_records"]:
        check(item)
    output = plan_path.parent
    destination = output / "inference"
    if destination.exists() or (output / "execution.json").exists():
        raise FileExistsError("Existing inference or execution report; no automatic resume")
    destination.mkdir()
    env = dict(os.environ, **plan["environment_overrides"])
    for key in plan["environment_unset"]:
        env.pop(key, None)
    argv = ["/usr/bin/time", "-v", "-o", str(output / "time.txt"), *plan["command"]]
    result = dict(status="running", plan=record(plan_path), command=argv,
                  job_id=os.environ["SLURM_JOB_ID"], attempts=1, scientific_scores_admitted=False)
    try:
        with (output / "inference.log").open("x") as stream:
            child = subprocess.Popen(argv, cwd=output, env=env, stdout=stream,
                                     stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = child.wait(timeout=23 * 3600)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                raise
        result["returncode"] = code
        if code:
            raise RuntimeError(f"Installed pipeline exited {code}")
        for item in plan["checked_records"]:
            check(item)
        result["status"] = "native_completed_pending_independent_scientific_readback"
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        result["logs"] = [record(output / name) for name in ("inference.log", "time.txt") if (output / name).exists()]
        with (output / "execution.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--prepare", type=Path)
    parser.add_argument("--run", type=Path)
    args = parser.parse_args()
    if args.prepare and args.repo and not args.run:
        prepare(args.repo.resolve(), args.prepare.absolute())
    elif args.run and not args.prepare:
        run(args.run.resolve())
    else:
        parser.error("Choose --repo ROOT --prepare OUTPUT or --run PLAN")
