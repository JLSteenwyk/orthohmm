"""Verify recovered BPO evidence for a future native inference launcher."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_blast_recovery_bpo import verify_admission
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_blast import PLAN_SHA
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMITTER_SHA = "b1b21e5e0a332c31bdf9b1b21341704a411ad1ecf81b44fd00b32dcb3122c934"


def completed(job):
    accounting = subprocess.check_output([
        "sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    row = require_completed_job(accounting, job)
    if any(row.get(k) != v for k, v in dict(
            NodeList="bizon", AllocCPUS="2", ReqMem="64G").items()):
        raise ValueError("Require completed two-CPU 64-GiB recovery validation")
    return row, accounting


def admitted_inputs(admission, root, species_record):
    if (admission["status"] != "recovered_orthomcl_bpo_checkpoint_admitted"
            or admission["checkpoint_admitted"] is not True
            or any(admission[k] is not False for k in (
                "accuracy_admitted", "publication_ready", "downstream_execution_authorized"))):
        raise ValueError("Require independently admitted recovered BPO checkpoint")
    content = admission["validation"]["content"]
    coverage = admission["query_coverage"]
    if content["input_proteins"] != 984137 or coverage["input_proteins"] != 984137:
        raise ValueError("Wrong recovered input universe")
    for actual, expected in (("source_hsp_rows", "hsp_rows"),
                             ("source_pair_blocks", "distinct_directed_pairs")):
        if content[actual] != coverage[expected]:
            raise ValueError("Recovered BPO counts differ from search coverage")
    base = root / "benchmarks/results/qfo_blast_recovery_bpo_v1/checkpoint"
    paths = dict(bpo=base / "all.bpo", offsets=base / "indexes/all_bpo.idx",
                 ranges=base / "indexes/all_bpo.se",
                 species=root / "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.gg")
    native = admission["native_inputs"]
    if (len(native) != 3 or {r["path"] for r in native}
            != {str(paths[k]) for k in ("bpo", "offsets", "ranges")}):
        raise ValueError("Wrong recovered native input inventory")
    inputs = {}
    for key, path in paths.items():
        matches = [r for r in admission["checked_records"] if r["path"] == str(path)]
        if key == "species":
            if species_record["path"] != str(path):
                raise ValueError("Wrong prepared species mapping path")
            matches.append(species_record)
        if (not matches or any(r != matches[0] for r in matches)
                or type(matches[0]["bytes"]) is not int or matches[0]["bytes"] <= 0):
            raise ValueError("Missing/conflicting admitted native input: " + key)
        if key != "species" and matches[0] not in native:
            raise ValueError("Native artifact differs from checked record")
        inputs[key] = matches[0]
    return inputs


def verify_inputs(root, path, digest, job, executor, commit):
    expected = root / "benchmarks/results/qfo_blast_recovery_bpo_admission_v1/report.json"
    if path != expected:
        raise ValueError("Unexpected recovered BPO admission path")
    if not executor.resolve().is_relative_to(root / "benchmarks/work"):
        raise ValueError("Require retained recovery admission executor")
    scheduler, accounting = completed(job)
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Changed recovered BPO admission executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--",
                    "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/admit_blast_recovery_bpo.py")
    if source["sha256"] != ADMITTER_SHA:
        raise ValueError("Unreviewed recovered BPO admission source")
    admission = read_frozen(path, digest)
    if admission["source"] != source:
        raise ValueError("Recovered BPO report differs from validator source")
    plan_path = root / "benchmark_tools/results/qfo_corrected_orthomcl_prepared_20260918.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    species_path = root / "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.gg"
    species = [r for r in plan["prepared_inputs"] if r["path"] == str(species_path)]
    if len(species) != 1:
        raise ValueError("Require unique species mapping from frozen input preparation")
    inputs = admitted_inputs(admission, root, species[0])
    parent_scheduler, _ = completed(admission["scheduler"]["JobIDRaw"])
    if parent_scheduler != admission["scheduler"]:
        raise ValueError("Recovered preparation accounting changed")
    search_record = admission["recovered_search"]
    search, _, search_checked, _, _ = verify_admission(
        root, Path(search_record["path"]), search_record["sha256"])
    if (search["query_coverage"] != admission["query_coverage"]
            or search_record not in admission["checked_records"]):
        raise ValueError("Recovered search identity or failure coverage changed")
    checked = [record(path), source, record(plan_path), *admission["checked_records"],
               *admission["validation"]["outputs"], *inputs.values(), *search_checked]
    for item in checked:
        check(item)
    return dict(status="recovered_native_input_evidence_verified_no_execution",
                admission=record(path), admission_job=str(job), scheduler=scheduler,
                accounting=accounting, executor=str(executor), commit=commit,
                inputs=inputs, index_validation=admission["validation"]["index_validation"],
                query_coverage=admission["query_coverage"], checked_records=checked,
                execution_authorized=False, accuracy_admitted=False, publication_ready=False,
                limitations=[
                    "Evidence gate only; native launcher must verify runtimes and stage fresh copies.",
                    "Admission job is supplied explicitly and checked in Slurm; report has no embedded admission job ID.",
                    "Search failures remain retained; successful conversion does not restore missing hits."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "admission", "executor"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    result = verify_inputs(args.root.resolve(), args.admission.resolve(),
                           args.admission_sha256, args.job, args.executor.resolve(), args.commit)
    print(json.dumps({"status": result["status"], "execution_authorized": False}))
