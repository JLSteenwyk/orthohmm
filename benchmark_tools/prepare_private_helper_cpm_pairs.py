"""Losslessly convert admitted recovered private high-CPM native pairs for QfO."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import admit_private_helper_cpm_phylogeny as native_checker
from benchmark_tools.admit_qfo_corrected_factorial_cell import gene_ownership
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_qfo_parameter_pairs import convert
from benchmark_tools.run_simulation_methods import read_frozen

PYTHON = "/home/bizon/anaconda3/bin/python"
ADMITTER_COMMIT = "d3c4c5c3525fcbeea97fa6b6d52cce05a0d5de29"
ADMITTER_SHA = "a21036cf1b916b65053c9363d19ed86bcc596a89fe2239a72da0ad5a56ea2666"
ADMITTER_EXECUTOR = "benchmarks/work/qfo_private_cpm_native_admission_executor_20261001"
ADMISSION = "benchmarks/work/qfo_private_cpm_native_admission_20261001.json"
NATIVE_PROTOCOL_SHA = "2e6db7249a8dddebe02674c955e7654457325eb756c7aeaf77d387b4e1c74542"
NATIVE_SUBMISSION_SHA = "3d3dde8a7f29d78a640655e6ac84abc44d9e2f89311327133b82ffe0fb5a70cd"
PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_CPM_PAIR_PROTOCOL_20261001.md"
OUTPUT = "benchmarks/results/qfo_cpm_private_pairs_v1/cpm_high"


def completed_admission(text, job):
    if type(job) is not str or not job.isascii() or not job.isdecimal() or job.startswith("0"):
        raise ValueError("Require explicit standalone admission job identity")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    parent = [row for row in rows if row["JobID"] == job]
    if len(parent) != 1 or tuple(parent[0][key] for key in
            ("JobIDRaw", "State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem")) != (
                job, "COMPLETED", "0:0", "bizon", "2", "64G"):
        raise ValueError("Require successfully completed private native admission")
    steps = [row for row in rows if row["JobID"].startswith(job + ".")]
    if (not any(row["JobID"] == job + ".batch" for row in steps)
            or len({row["JobID"] for row in steps}) != len(steps)
            or any(row["State"] != "COMPLETED" or row["ExitCode"] != "0:0" for row in steps)):
        raise ValueError("Require successful terminal private native admission steps")
    return parent[0]


def accounting(job):
    return subprocess.check_output(["sacct", "-j", job, "-P",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)


def admission_submission(root, job, sha):
    path = root / f"benchmark_tools/results/qfo_private_cpm_native_admission_submission_{job}.json"
    submission = read_frozen(path, sha)
    executor = root / ADMITTER_EXECUTOR
    script = "benchmark_tools/results/qfo_private_cpm_native_admission_20261001.sh"
    if (submission["status"] != "private_recovered_qfo_high_cpm_native_admission_submitted"
            or submission["job_id"] != job or submission["executor"] != str(executor)
            or submission["executor_commit"] != ADMITTER_COMMIT or submission["executor_clean"] is not True
            or any(submission[key] is not False for key in
                   ("native_inference", "accuracy_evaluated", "controlled_timing", "publication_ready"))
            or submission["submission_argv"] != ["sbatch", "--parsable", str(executor / script),
                str(executor), ADMITTER_COMMIT, NATIVE_PROTOCOL_SHA]
            or submission["parent_submission"] != record(root / native_checker.SUBMISSION)
            or submission["parent_submission"]["sha256"] != NATIVE_SUBMISSION_SHA):
        raise ValueError("Private native admission submission identity differs")
    names = ("benchmark_tools/admit_private_helper_cpm_phylogeny.py", native_checker.PROTOCOL,
             script, "tests/unit/test_admit_private_helper_cpm_phylogeny.py")
    pins = [record(executor / name) for name in names]
    if submission["source_records"] != pins:
        raise ValueError("Private native admission submitted source inventory differs")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER_COMMIT:
        raise ValueError("Private native admission executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm",
                    "tests/unit/test_admit_private_helper_cpm_phylogeny.py"], check=True, capture_output=True)
    if (pins[0]["sha256"] != ADMITTER_SHA or pins[1]["sha256"] != NATIVE_PROTOCOL_SHA
            or record(native_checker.__file__)["sha256"] != ADMITTER_SHA):
        raise ValueError("Private native admission or reused producer verifier source changed")
    return pins[0], [record(path), submission["parent_submission"], *pins]


def validate_native(native, source, verified, cell, root):
    if (native["status"] != "private_recovered_cpm_native_pairs_verified_unscored"
            or native["arm"] != "cpm_high" or type(native["index"]) is not int or native["index"] != 1
            or native["source"] != source or native["cell"] != cell
            or native["candidate_admission"] != verified["candidates"]["admission_record"]
            or native["private_native_admission"] != verified["private_control"]["admission"]
            or native["protocol"] != record(root / native_checker.PROTOCOL)
            or native["protocol"]["sha256"] != NATIVE_PROTOCOL_SHA
            or native["submission"] != record(root / native_checker.SUBMISSION)
            or native["submission"]["sha256"] != NATIVE_SUBMISSION_SHA
            or any(native[key] is not False for key in
                   ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))
            or type(native["native_pair_count"]) is not int or native["native_pair_count"] <= 0
            or any(native[key] not in native["checked_records"] for key in ("native_pairs", "source", "protocol", "submission"))
            or native["native_group_integrity"]["native_manifest"] not in native["checked_records"]):
        raise ValueError("Wrong or incomplete recovered private native admission")


def prepare(root, admission_job, admission_sha, submission_sha, protocol_sha):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon" or os.environ.get("SLURM_MEM_PER_NODE") != "65536"
            or os.environ.get("SLURM_ARRAY_TASK_ID") or sys.executable != PYTHON):
        raise ValueError("Require scheduled standalone two-CPU/64-GiB bizon conversion controller")
    if Path.cwd().resolve() != root:
        raise ValueError("Run from original repository verification directory")
    output = root / OUTPUT
    if output.resolve() != output or output.exists() or output.is_symlink():
        raise FileExistsError("Private conversion output must be fresh and direct")
    observed = accounting(admission_job)
    scheduler = completed_admission(observed, admission_job)
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Private conversion protocol identity changed")
    source, submitted = admission_submission(root, admission_job, submission_sha)
    admission = record(root / ADMISSION)
    native = read_frozen(Path(admission["path"]), admission_sha)
    if admission["sha256"] != admission_sha:
        raise ValueError("Private native admission changed during read")
    producer, verified = native_checker.producer(root, root / native_checker.EXECUTOR)
    baseline, candidates = verified["baseline"], verified["candidates"]
    cell, _, _, _ = producer.native_command(verified, root / producer.OUTPUT)
    validate_native(native, source, verified, cell, root)
    environment_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(environment_path, ENV_SHA)
    mappings = [ref for ref in environment["reference_files"] if Path(ref["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen QfO reference mapping")
    mapping = mappings[0]
    helpers = [record(module.__file__) for name, module in sorted(sys.modules.items())
        if name.startswith(("benchmark_tools.", "qfo_benchmark.")) and getattr(module, "__file__", None)]
    checked = [record(__file__), protocol, *submitted, admission, *native["checked_records"],
               *verified["checked_records"], record(environment_path), mapping, *helpers]
    for ref in checked:
        check(ref)
    output.mkdir(parents=True, exist_ok=False)
    recheck = output / "native_admission_recheck.json"
    command = [PYTHON, "-B", source["path"], "--root", str(root),
        "--submission-sha256", NATIVE_SUBMISSION_SHA, "--protocol-sha256", NATIVE_PROTOCOL_SHA,
        "--output", str(recheck)]
    report = dict(status="preparing", arm="cpm_high", index=1, participant="ohmm_qfo_parameter_cpm_high",
        job_id=os.environ["SLURM_JOB_ID"], admission_job=admission_job, admission_scheduler=scheduler,
        admission_accounting=observed, native_admission=admission, native_input=native["native_pairs"],
        candidate_admission=candidates["admission_record"], private_native_admission=verified["private_control"]["admission"],
        protocol=protocol, source=record(__file__), helpers=helpers, native_recheck_command=command,
        mapping=mapping, input_fastas=baseline["manifest"]["input_fastas"], checked_records=checked,
        semantics="native phylogenetically inferred pairs", accuracy_evaluated=False, scoring_admitted=False,
        controlled_timing=False, publication_ready=False)
    with (output / "preflight.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    checked.append(record(output / "preflight.json"))
    try:
        with (output / "native_recheck.log").open("x") as log:
            subprocess.run(command, check=True, stdout=log, stderr=log, cwd=root)
        fresh = record(recheck)
        if read_frozen(recheck, fresh["sha256"]) != native:
            raise ValueError("Fresh private native admission differs from retained report")
        checked.extend([fresh, record(output / "native_recheck.log")])
        metadata = native["native_group_integrity"]["native_manifest"]
        owners, _ = gene_ownership(baseline["manifest"], read_frozen(Path(metadata["path"]), metadata["sha256"]),
                                  Path(candidates["arm"]["partition"]["path"]))
        counts = convert(native, owners, Path(mapping["path"]), output)
        partial = [record(output / name) for name in ("pairs.partial.tsv", "pairs.qfo.partial.tsv")]
        if any(partial[0][key] != partial[1][key] for key in ("bytes", "sha256")):
            raise ValueError("Private QfO filtering changed converted pair bytes")
        checked.extend(partial)
        for ref in checked:
            check(ref)
        if producer.verify_sources(root) != verified or completed_admission(accounting(admission_job), admission_job) != scheduler:
            raise ValueError("Private conversion evidence/runtime/admission completion changed")
        for ref in checked:
            check(ref)
        for ref, name in zip(partial, ("pairs.tsv", "pairs.qfo.tsv")):
            Path(ref["path"]).rename(output / name)
            checked.remove(ref)
        final = [record(output / name) for name in ("pairs.tsv", "pairs.qfo.tsv")]
        if any(any(before[key] != after[key] for key in ("bytes", "sha256")) for before, after in zip(partial, final)):
            raise ValueError("Private converted pair bytes changed during publication")
        checked.extend([*final, record(output / "conversion_counts.json")])
        report.update(status="private_recovered_cpm_native_pairs_prepared_unscored", **counts,
            pairs=final[0], filtered_pairs=final[1], native_admission_recheck=fresh,
            conversion_counts=record(output / "conversion_counts.json"))
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--admission-job", required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-submission-sha256", required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.admission_job, args.admission_sha256,
            args.admission_submission_sha256, args.protocol_sha256)
