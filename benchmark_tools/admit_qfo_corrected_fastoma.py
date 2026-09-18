"""Independently admit corrected FastOMA native evidence before pair scoring."""

import argparse
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_fastoma_orthoxml import read_xml, check_root_table, check_pair_scope
from benchmark_tools.audit_fastoma_tasks import audit_tasks, option
from benchmark_tools.fastoma_to_pairwise import input_owners
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_fastoma import preflight
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "7aa174c3f6bf244144c1940bfa32db2e3ca80d7b"


def validate_execution(execution, expected, scheduler, runner):
    if execution.get("status") != "process_succeeded_pending_native_admission" or execution.get("exit_code") != 0:
        raise ValueError("FastOMA execution has not succeeded")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "180"
            or scheduler["ReqMem"] != "720G" or execution["job_id"] != scheduler["JobIDRaw"]
            or execution["node"] != "bizon" or execution["source"] != runner):
        raise ValueError("Wrong FastOMA source, scheduler identity or allocation")
    if any(execution[key] != value for key, value in expected.items()):
        raise ValueError("Execution differs from freshly verified preflight")
    if any(execution[key] is not False for key in ("accuracy_admitted", "native_outputs_validated", "publication_ready")):
        raise ValueError("Unexpected execution admission state")
    start, end = execution["started_epoch"], execution["finished_epoch"]
    if not math.isfinite(start) or not math.isfinite(end) or not end >= start > 0:
        raise ValueError("Invalid native execution timestamps")
    output = Path(expected["output_root"])
    for key, path in (("log", output / "native.log"), ("timing", output / "time.txt"),
                      ("trace", output / "run/trace.txt"), ("nextflow_log", output / "run/.nextflow.log")):
        if execution[key]["path"] != str(path) or execution[key]["bytes"] <= 0:
            raise ValueError("Unexpected or empty execution artifact")


def tree_clades(path, species):
    tree = Phylo.read(path, "newick")
    leaves = [node.name for node in tree.get_terminals()]
    if len(leaves) != len(set(leaves)) or set(leaves) != set(species):
        raise ValueError("Checked FastOMA tree species differ")
    for node in tree.find_clades():
        if node.branch_length is not None and not math.isfinite(node.branch_length):
            raise ValueError("Nonfinite species tree branch")
    return {frozenset(leaf.name for leaf in node.get_terminals()) for node in tree.get_nonterminals()}


def bind_task_outputs(tasks, output, inputs):
    selected = {t["trace"]["name"].split(" (", 1)[0]: Path(t["directory"]) for t in tasks["tasks"]}
    bindings = []
    names = {"collect_subhogs": ("FastOMA_HOGs.orthoxml", "RootHOGs.tsv", "OrthologousGroups.tsv"),
             "extract_pairwise_ortholog_relations": ("orthologs.tsv.gz",),
             "check_input": ("species_tree_checked.nwk",), "fastoma_report": ("report.ipynb", "report.html")}
    for process, files in names.items():
        for name in files:
            native, published = record(selected[process] / name), record(output / name)
            if any(native[key] != published[key] for key in ("bytes", "sha256")):
                raise ValueError("Published artifact differs from native task: " + name)
            bindings.append({"native": native, "published": published})
    pair_input = record(selected["extract_pairwise_ortholog_relations"] / "FastOMA_HOGs.orthoxml")
    published = record(output / "FastOMA_HOGs.orthoxml")
    if any(pair_input[key] != published[key] for key in ("bytes", "sha256")):
        raise ValueError("Pair extraction used different OrthoXML")
    bindings.append({"native": pair_input, "published": published})
    proteomes = {Path(row["path"]).name: row for row in inputs}
    for task in tasks["tasks"]:
        if task["trace"]["name"].split(" (", 1)[0] != "omamer_run":
            continue
        directory = Path(task["directory"])
        query = option(shlex.split((directory / ".command.sh").read_text(), comments=True), "--query")
        staged = proteomes[query]
        native = record(directory / query)
        if any(native[k] != staged[k] for k in ("bytes", "sha256")):
            raise ValueError("OMAmer used a different proteome")
        bindings.append({"native": native, "staged": staged})
        mapping = record(directory / (query + ".hogmap"))
        for path in (output / "hogmap" / (query + ".hogmap"),
                     selected["infer_roothogs"] / "hogmaps" / (query + ".hogmap")):
            observed = record(path)
            if any(mapping[k] != observed[k] for k in ("bytes", "sha256")):
                raise ValueError("OMAmer mapping changed at publication or root inference")
            bindings.append({"native": mapping, "consumed_or_published": observed})
    return bindings


def admit(root, job, stage_job, destination):
    if destination.exists():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    executor = root / "benchmarks/work/publication_qfo_corrected_fastoma_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Inference executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    run_root = root / "benchmarks/results/qfo_corrected_fastoma_v1"
    execution_path = run_root / "execution.json"
    execution_record = record(execution_path)
    execution = json.loads(execution_path.read_text())
    stage_path = root / "benchmarks/work/qfo_corrected_fastoma_staging_20260918.json"
    if execution["stage"]["path"] != str(stage_path):
        raise ValueError("Unexpected staging manifest path")
    expected, _ = preflight(root, stage_path, execution["stage"]["sha256"], stage_job)
    validate_execution(execution, expected, scheduler, record(executor / "benchmark_tools/run_qfo_corrected_fastoma.py"))
    stage = read_frozen(stage_path, execution["stage"]["sha256"])
    output = run_root / "output"
    observed = [record(p) for p in sorted(output.rglob("*")) if p.is_file()]
    if not observed or observed != execution["outputs"]:
        raise ValueError("Published native output inventory changed")
    inputs = [row["staged"] for row in stage["copies"] if Path(row["staged"]["path"]).suffix == ".fa"]
    tasks = audit_tasks(run_root / "run/trace.txt", run_root / "work", [Path(r["path"]).name for r in inputs])
    bindings = bind_task_outputs(tasks, output, inputs)
    owners = input_owners([Path(r["path"]) for r in inputs])
    if len(owners) != 984137 or len(set(owners.values())) != 78:
        raise ValueError("Corrected input universe differs")
    tree = output / "species_tree_checked.nwk"
    if tree_clades(tree, owners.values()) != tree_clades(Path(stage["input_directory"]) / "species_tree.nwk", owners.values()):
        raise ValueError("FastOMA changed supplied species-tree topology")
    content, members = read_xml(output / "FastOMA_HOGs.orthoxml", owners)
    content["root_table_members"] = check_root_table(output / "RootHOGs.tsv", members)
    content.update(check_pair_scope(output / "orthologs.tsv.gz", owners, members))
    checked = [execution_record, execution["source"], *expected["checked_records"], *observed,
               *[execution[key] for key in ("log", "timing", "trace", "nextflow_log")],
               *[r for task in tasks["tasks"] for r in task["files"]],
               *[r for binding in bindings for r in binding.values()]]
    for item in checked:
        check(item)
    if {str(p.resolve()) for p in output.rglob("*") if p.is_file()} != {r["path"] for r in observed}:
        raise ValueError("Output membership changed during admission")
    report = {"status": "corrected_fastoma_native_evidence_admitted", "source": record(__file__),
              "scheduler": scheduler, "accounting": accounting, "execution": execution_record,
              "checked_records": checked, "tasks": tasks, "task_output_bindings": bindings,
              "content": content, "native_pairs": record(output / "orthologs.tsv.gz"),
              "input_fastas": inputs, "accuracy_evaluated": False, "publication_ready": False,
              "helpers": [record(Path(__file__).with_name(name)) for name in (
                  "run_qfo_corrected_fastoma.py", "audit_fastoma_tasks.py", "audit_fastoma_orthoxml.py", "fastoma_to_pairwise.py")],
              "limitations": ["Supplied corrected OrthoFinder tree, not independent FastOMA tree inference.",
                  "Native integrity and scope checks, not biological accuracy or reference-relative coverage.",
                  "Failed/retried task chains require separate review; no automatic exception applied.",
                  "Distinct native pair conversion and independent QfO scoring still required.",
                  "Shared-host inference is not matched dedicated timing."]}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--stage-job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.stage_job, args.output.resolve())
