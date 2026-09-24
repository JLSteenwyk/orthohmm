"""Use the admitted recovered seed with unchanged candidate expansion settings."""

import argparse
import importlib
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recovered_cpm_seed_evidence import evidence

BUILD_SHA = "b6f3325e2e6c33eecbf74256e2bd88d339015188e47c9296fedea1155029b978"


def prepare(root):
    from benchmark_tools.prepare_qfo_candidate_neighborhood import require_environment
    require_environment(os.environ)
    if os.environ.get("SLURM_MEM_PER_NODE") != "65536":
        raise ValueError("Require 64-GiB recovery candidate allocation")
    output = root / "benchmarks/results/qfo_cpm_recovered_candidates_v1"
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    admitted = evidence(root)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    from benchmark_tools.prepare_qfo_cpm_candidates import build, BASELINE_SHA
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.cpm_replay_context import evidence as context_evidence, REPLAY_SHA
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.run_blast_recovery_batch import save_status

    builder = record(Path(__file__).with_name("prepare_qfo_cpm_candidates.py"))
    if builder["sha256"] != BUILD_SHA:
        raise ValueError("Candidate expansion wrapper changed")
    plan, plan_record, _, _ = corrected_evidence(root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json", REPLAY_SHA)
    context = context_evidence(root, plan, plan_record, "cpm_high")
    baseline_path = root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    parameters = baseline["candidate_arms"]["p1_c1"]["expansion"]["parameters"]
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == baseline["runtime_before"] == baseline["runtime_after"]:
        raise ValueError("Recovered candidate runtime differs")
    sys.path.insert(0, str(launcher))
    try:
        engine = importlib.import_module("orthohmm.orthohmm")
        accuracy = importlib.import_module("orthohmm.accuracy")
    finally:
        sys.path.remove(str(launcher))
    if any(Path(module.__file__).resolve().parent != launcher / "orthohmm" for module in (engine, accuracy)):
        raise ValueError("Wrong frozen candidate scientific import")
    # The numeric auditor imports accuracy, so pin the scientific package first.
    from benchmark_tools.audit_accuracy_checkpoint import audit
    from benchmark_tools.audit_candidate_arm import audit as audit_arm
    checkpoint_record = plan["checkpoint_manifest"]
    checkpoint = Path(checkpoint_record["path"]).parent
    numeric = audit(checkpoint, checkpoint_record["sha256"])
    if numeric["summary"] != baseline["numeric_checkpoint"]["summary"] or numeric["summary"]["species"] != 78:
        raise ValueError("Changed candidate numeric checkpoint")
    names, species, queries, targets, scores = accuracy.load_accuracy_checkpoint(checkpoint, verify=False)
    if len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong recovered candidate gene universe")
    helpers = [record(path) for path in sorted(Path(__file__).parent.glob("*.py"))]
    records = [plan_record, record(baseline_path), record(runtime_path), builder,
               *context["checked_records"], *admitted["checked_records"], *helpers]
    for item in records:
        check(item)
    output.mkdir()
    report = dict(status="recovered_candidates_preparing", arm="cpm_high", context=context,
        source=record(__file__), builder=builder, recovery_admission=admitted, inputs=records,
        runtime_before=runtime, numeric_checkpoint=numeric, job_id=os.environ["SLURM_JOB_ID"],
        executor_commit=subprocess.check_output(["git", "-C", str(Path(__file__).resolve().parent.parent),
                                                "rev-parse", "HEAD"], text=True).strip(),
        accuracy_evaluated=False, publication_ready=False)
    save_status(output / "manifest.json", report)
    try:
        report.update(build(engine, parameters, admitted["seed_partition"], names, species,
                            (queries, targets, scores), output / "candidate", audit_arm))
        if audit(checkpoint, checkpoint_record["sha256"]) != numeric:
            raise ValueError("Checkpoint changed during recovered candidate construction")
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Runtime changed during recovered candidate construction")
        for item in [*records, *report["candidate_arm"]["output_files"]]:
            check(item)
        report.update(status="recovered_cpm_candidates_prepared_pending_admission", limitations=[
            "Candidate logic and thresholds unchanged; seed is independently admitted checkpoint recovery.",
            "Independent candidate admission and downstream phylogeny/scoring remain required.",
            "Shared-host incremental construction is not controlled end-to-end timing."])
    except BaseException as error:
        report.update(status="recovered_candidate_preparation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "manifest.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve())
