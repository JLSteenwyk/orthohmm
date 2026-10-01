"""Isolated stdlib readback of the completed, unscored private QfO admission."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess

COMMIT = "9f0fcb68b849aa8af6180c73c9ba82d3e4f4b8b1"
ADMISSION = "benchmarks/work/qfo_private_phylogeny_control_admission_v2_20261001.json"
ADMISSION_SHA = "6865539281f9cb8c89c04b8ca16e2277056b4cc2b3293c330a12a3115879af97"
NATIVE_FILES = {
    "orthohmm_root_hogs.tsv", "orthohmm_pairwise_orthologs.tsv",
    "orthohmm_pairwise_orthologs_confidence.tsv", "orthohmm_reconciliation_nodes.tsv",
    "orthohmm_hierarchical_orthogroups.tsv", "species_tree.rooted.nwk",
}


def record(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def completed(text):
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    parents = []
    for job, cpus, memory in (("22387", "32", "192G"), ("22389", "2", "64G")):
        selected = [r for r in rows if r["JobID"] == job]
        if len(selected) != 1 or tuple(selected[0][k] for k in
                ("JobIDRaw", "State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
                    job, "COMPLETED", "0:0", cpus, memory, "bizon"):
            raise ValueError("Require exact completed private control/admission parents")
        steps = [r for r in rows if r["JobID"].startswith(job + ".")]
        if not any(r["JobID"] == job + ".batch" for r in steps) or any(
                r["State"] != "COMPLETED" or r["ExitCode"] != "0:0" for r in steps):
            raise ValueError("Require successful terminal private control/admission steps")
        parents.append(selected[0])
    return parents


def readback(root):
    command = ["sacct", "-j", "22387,22389", "-P",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"]
    accounting = subprocess.check_output(command, text=True)
    scheduler = completed(accounting)
    pin = record(root / ADMISSION)
    if pin["sha256"] != ADMISSION_SHA:
        raise ValueError("Private admission report changed")
    report = json.loads(Path(pin["path"]).read_bytes())
    if (report["status"] != "private_qfo_phylogeny_deployment_admitted_unscored"
            or report["job_id"] != "22387" or report["scheduler"] != scheduler[0]
            or report["recovered_cpm_inference_authorized"] is not True
            or any(report[k] is not False for k in
                ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))
            or type(report["native_pair_count"]) is not int or report["native_pair_count"] != 5959560
            or set(report["native_comparison"]) != NATIVE_FILES):
        raise ValueError("Private admission has unexpected semantics or native outputs")
    submission_pin = record(root / "benchmark_tools/results/qfo_private_phylogeny_admission_submission_22389.json")
    submission = json.loads(Path(submission_pin["path"]).read_bytes())
    executor = root / "benchmarks/work/qfo_private_phylogeny_admission_v2_executor_20261001"
    if (submission["job_id"] != "22389" or submission["executor"] != str(executor)
            or submission["executor_commit"] != COMMIT or submission["executor_clean"] is not True
            or report["source"] != submission["source_records"][0]
            or subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT):
        raise ValueError("Wrong private admission source/submission")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"],
                   check=True, capture_output=True)
    bindings = []
    for item in submission["source_records"]:
        path = Path(item["path"])
        blob = subprocess.check_output(["git", "-C", str(executor), "show", f"{COMMIT}:{path.relative_to(executor)}"])
        if record(path) != item or blob != path.read_bytes():
            raise ValueError("Private admission source/Git identity changed")
        bindings.append({**item, "git_revision": COMMIT})
    records = [pin, record(__file__), submission_pin, submission["parent_submission"],
        submission["preserved_failed_attempt"], *submission["source_records"], *report["checked_records"],
        record(root / "benchmarks/work/qfo_private_phylogeny_admission_22389.log")]
    for item in report["native_comparison"].values():
        old, new = item["historical"], item["current"]
        if item["byte_equal"] is not True or (old["bytes"], old["sha256"]) != (new["bytes"], new["sha256"]):
            raise ValueError("Private native output parity disagrees")
        records.extend([old, new])
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting readback file identity")
        unique[item["path"]] = item
    for _ in range(2):
        for item in unique.values():
            if record(item["path"]) != item:
                raise ValueError("Readback bound file changed: " + item["path"])
    if completed(subprocess.check_output(command, text=True)) != scheduler:
        raise ValueError("Private completion changed during readback")
    return dict(status="private_qfo_phylogeny_admission_independently_read_back_unscored",
        source=record(__file__), admission=pin, submission=submission_pin, scheduler=scheduler,
        accounting=accounting, source_revision=COMMIT, source_git_bindings=bindings,
        checked_records_reverified=len(unique), native_pair_count=report["native_pair_count"],
        cache_use=report["cache_use"], native_comparison=report["native_comparison"],
        recovered_cpm_inference_authorized=True, accuracy_evaluated=False, scoring_admitted=False,
        controlled_timing=False, publication_ready=False,
        limitations=["Isolated stdlib file/completion/source readback, not a third native partition/pair/tree validation.",
            "Complete native inventory and validators are bound through the unchanged completed admission.",
            "Native byte identity under checkpoint reuse, not fresh all-tree inference, reference accuracy or controlled timing."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(readback(args.root.resolve()), indent=2, sort_keys=True, allow_nan=False))
