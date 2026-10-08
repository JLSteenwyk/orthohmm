"""Compose original scientific/resource kernels with explicit late-runtime review."""

import argparse
import json
import os
from pathlib import Path
import pwd
import re
import subprocess
import time

from benchmark_tools import review_native11_postterminal_runtime as runtime_kernel
from benchmark_tools import review_native11_postterminal_runtime_v2 as runtime_v2
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.native_factorial_allocated_execution import (
    amendment, validate_request, verify_terminal, native_command, MEMORY, SCOPE)
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_allocated_threadripper_scaling import replay
from benchmark_tools.review_allocated_native_factorial_attempt import bind_session, environment_review
from benchmark_tools.review_native_factorial_attempt import save, scheduler_fields, classify, resource_review, DISCLOSURE
from benchmark_tools.run_native11_fault_reported_review import diagnostic_gate
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import Evidence, require
from benchmark_tools.verify_lineage_native_provenance import same


SCHEMA = "native11_composed_terminal_review_v1"
DESTINATION = ROOT / "benchmarks/work/native11_composed_terminal_review_20261008_v1"
BATCH = ROOT / "benchmark_tools/results/native11_composed_terminal_review_20261008_v1.sh"
RUNTIME = ROOT / "benchmarks/work/native11_postterminal_runtime_review_20261008_v2/runtime.json"
RUNTIME_SHA = "ed8ca2d27e3f85d3fdf14952d33b2db7d64dca0bb76d5113fb483bfe7e27a8c6"
OUTPUT = ROOT / "benchmarks/work/native11_standalone_diagnostic_20261008_v1/outputs.json"
OUTPUT_SHA = "2d48c2f447c72d04e3e96c9094fe4a941404bd19b8daff56aa550ac3fee712be"
READBACK = ROOT / "benchmark_tools/results/native11_standalone_diagnostic_terminal_24031_20261008_v1.json"
READBACK_SHA = "e0aa4821e0b6295be7830ddf171e6de33921efd47d2b9759631054e2137a0f5a"
EXTRA_SOURCES = {
    "benchmark_tools/run_native11_fault_reported_review.py":
        "2e7f9928de9e5f417bf7d842468fb827a8779f18318f426a21e9664474f3e6cb",
    "benchmark_tools/review_native11_postterminal_runtime.py": runtime_v2.KERNEL_SHA,
    "benchmark_tools/review_native11_postterminal_runtime_v2.py":
        "46f6bc315fc41dbfd08979dd8b8ee7514344ad19a58a15af6c30c70919116e54",
}


def allocation_gate(raw, job, digest):
    lines = [line for line in raw.splitlines() if line.strip()]
    require(len(lines) == 1, "Require one composed-review controller record")
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", lines[0])
    fields = dict(pairs)
    require(len(fields) == len(pairs), "Duplicate controller field")
    expected = dict(JobId=str(job), JobName="ohmm_native11_composed", JobState="RUNNING",
        Partition="gpu", NodeList="bizon", NumNodes="1", NumCPUs="2", NumTasks="1",
        MinMemoryNode="128G", TimeLimit="06:00:00", Requeue="0", Restarts="0",
        Command=str(BATCH), WorkDir=str(ROOT), Comment=digest,
        UserId=f"{pwd.getpwuid(os.getuid()).pw_name}({os.getuid()})")
    expected["CPUs/Task"] = "2"
    require(type(job) is int and job > 24033
        and all(fields.get(k) == value for k, value in expected.items())
        and not {"ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset"}.intersection(fields),
        "Composed review allocation/owner/source binding differs")
    return fields


def runtime_gate(report, request_ref, request):
    require(report.get("schema") == runtime_v2.SCHEMA
        and report.get("status") == "historical_brackets_and_current_exact_additions_revalidated"
        and report.get("source") == record(runtime_v2.__file__)
        and report.get("runtime_kernel_source") == record(runtime_kernel.__file__)
        and report.get("request") == request_ref and report.get("plan") == request["plan"]
        and report.get("amendment") == request["amendment"]
        and report.get("job_id") == 23985 and report.get("index") == 11
        and report.get("cell") == "p1_c1_r0"
        and all(report.get(k) is True for k in ("historical_first_terminal_rechecked",
            "current_original_entries_unchanged", "current_private_inventory_equality"))
        and all(report.get(k) is False for k in ("current_original_inventory_equality",
            "continuous_runtime_integrity_established", "full_review_admitted", "terminal_reviewed",
            "next_identity_authorized", "accuracy_evaluated", "automatic_retry", "publication_ready")),
        "Require truthfully labeled successful postterminal runtime component")


def output_binding(request_ref, request, execution, run, evidence):
    output_ref, readback_ref = record(OUTPUT), record(READBACK)
    require(output_ref["sha256"] == OUTPUT_SHA and readback_ref["sha256"] == READBACK_SHA,
        "Standalone semantic component/readback changed")
    output, readback = read(output_ref), read(readback_ref)
    require(readback.get("status") == "standalone_semantic_diagnostic_completed_bound"
        and readback.get("output") == output_ref and readback.get("job_id") == 24031
        and readback.get("standalone_semantics_passed") is True
        and readback.get("full_review_admitted") is False, "Standalone independent readback differs")
    raw, fields = accounting(24031, include_memory=True)
    diagnostic_gate(output, fields, request_ref, execution["historical_plan"], request["amendment"],
        record(ROOT / "benchmark_tools/validate_allocated_native_factorial_outputs.py"))
    require(output.get("semantic_validator_source") == record(ROOT / "benchmark_tools/validate_native_factorial_outputs.py")
        and output.get("execution_scope") == SCOPE
        and output.get("allocated_ready") == record(Path(run["output_root"]) / "measurement/ready.json"),
        "Original semantic source or placement binding differs")
    refs = [output_ref, readback_ref, readback["submission"], readback["release"],
        *readback["direct_invocation_files"].values(), *readback["independent_checked_records"],
        *output["checked_files"], *output["evidence"], output["source"], output["semantic_validator_source"]]
    for ref in refs:
        check(ref)
        evidence.bind(ref["path"])
    return output_ref, output, dict(producer_job_id=24031, raw_accounting=raw, fields=fields,
        semantic_validator_reexecuted=False, readback=readback_ref)


def review(expected_source_sha256):
    source = record(__file__)
    require(source["sha256"] == expected_source_sha256, "Prospective composed reviewer changed")
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2"
        and os.environ.get("SLURM_JOB_ID", "").isdigit(), "Require separate scheduled two-CPU review")
    job = int(os.environ["SLURM_JOB_ID"])
    require(not DESTINATION.exists() and DESTINATION.resolve() == DESTINATION,
        "Require fresh direct composed-review destination")
    DESTINATION.mkdir(parents=True, exist_ok=False)
    evidence = Evidence()
    try:
        argv = ["scontrol", "show", "job", str(job), "--oneliner"]
        observation = subprocess.run(argv, capture_output=True, text=True, check=True, timeout=10)
        allocation = allocation_gate(observation.stdout, job, expected_source_sha256)
        allocation_ref = save(DESTINATION / "allocation.json", dict(command=argv,
            stdout=observation.stdout, stderr=observation.stderr, fields=allocation))
        request_ref = record(runtime_kernel.REQUEST)
        require(request_ref["sha256"] == runtime_kernel.REQUEST_SHA, "Native request changed")
        request = read(request_ref)
        execution, plan = amendment(request["amendment"])
        validate_request(request, request["amendment"], execution, 23985)
        run = plan["runs"][11]
        require(request["index"] == 11 and run["cell"] == "p1_c1_r0", "Wrong native identity")
        for path, digest in {**runtime_kernel.SOURCES, **EXTRA_SOURCES}.items():
            ref = record(ROOT / path)
            require(ref["sha256"] == digest, "Composed-review frozen kernel changed")
            evidence.bind(ref["path"])
        for ref in [source, record(BATCH), request_ref, request["amendment"], request["plan"],
                    plan["baseline"], *execution["new_sources"], *plan["helper_sources"], *plan["evidence"]]:
            check(ref)
            evidence.bind(ref["path"])
        terminal = verify_terminal(23985)
        runtime_v2.terminal_gate(terminal, request_ref)
        fields = scheduler_fields(terminal)
        scheduler_ref = save(DESTINATION / "scheduler.json", terminal)
        root = Path(run["output_root"])
        session = Path(plan["panel_root"]) / "sessions/run_11"
        result, verification = evidence.json(session / "result.json"), evidence.json(root / "verification.json")
        bind_session(result, verification, request_ref, request, request["amendment"], run)
        runtime_ref = record(RUNTIME)
        require(runtime_ref["sha256"] == RUNTIME_SHA, "Prospective runtime component changed")
        retained_runtime = read(runtime_ref)
        runtime_gate(retained_runtime, request_ref, request)
        for ref in [runtime_ref, *retained_runtime["evidence"]]:
            check(ref)
            evidence.bind(ref["path"])
        classification = read(retained_runtime["classification"])
        prior = read(retained_runtime["prior_runtime"])
        times = runtime_kernel.chronology(classification, read(classification["failed_review_observation"]))
        runtime = runtime_kernel.replay_runtime(plan, run, session, verification, prior, evidence,
            retained_runtime["current_additions"], times)
        runtime.update(schema="native11_composed_runtime_gate_v1", source=source,
            runtime_component=runtime_ref, runtime_kernel_source=record(runtime_kernel.__file__))
        fresh_runtime_ref = save(DESTINATION / "runtime.json", runtime)
        output_ref, outputs, producer = output_binding(request_ref, request, execution, run, evidence)
        producer_ref = save(DESTINATION / "semantic_producer.json", producer)
        baseline = read(plan["baseline"])
        replayed = replay(root / "measurement", 23985, native_command(request["amendment"], run, baseline))
        require(same(replayed["measured"], verification["measurement"]), "Collector and wrapper differ")
        require(classify(terminal, result, verification, replayed) == "native_success", "Native outcome differs")
        for ref in replayed["evidence"]:
            check(ref)
            evidence.bind(ref["path"])
        done = evidence.json(root / "measurement/done.json")
        resources = resource_review(replayed, done, 23985)
        resources_ref = save(DESTINATION / "resources.json", resources)
        environment = environment_review(request_ref, request, request["amendment"], run,
            result, root / "measurement", done, evidence)
        require(environment["sampled_environment_evidence_valid"] is True, "Shared-host evidence invalid")
        environment_ref = save(DESTINATION / "environment.json", environment)
        require(outputs["native_cpu_ids"] == replayed["native_cpu_ids"], "Output/resource CPU placement differs")
        # Retain the whole replay's raw evidence without serializing another multi-GB host-process expansion.
        replay_summary_ref = save(DESTINATION / "resource_replay_summary.json", dict(
            schema="native11_composed_replay_summary_v1", source=source,
            replay_source=replayed["source"], common_replay_source=replayed["common_replay_source"],
            full_original_replay_executed=True, native_outcome=replayed["native_outcome"],
            native_exit_code=replayed["native_exit_code"], evidence=replayed["evidence"],
            original_expanded_replay_serialized=False, scientific_timings_admitted=False))
        checked = evidence.finish()
        check(source)
        return save(DESTINATION / "review.json", dict(schema=SCHEMA, status="native_success",
            source=source, job_id=23985, review_job_id=job, index=11, cell=run["cell"],
            dataset=run["dataset"], repeat=run["repeat"], request=request_ref,
            plan=request["plan"], amendment=request["amendment"], scheduler=scheduler_ref,
            allocation=allocation_ref, scheduler_state=fields.get("JobState", fields.get("State")),
            scheduler_exit_code=fields["ExitCode"], evidence=checked,
            reviews=dict(runtime=fresh_runtime_ref, resources=resources_ref, environment=environment_ref,
                outputs_or_failure=output_ref), semantic_producer=producer_ref, resource_replay=replay_summary_ref,
            native_cpu_ids=replayed["native_cpu_ids"], allocated_placement=replayed["allocated_placement"],
            resources=resources["primary"], resource_scopes=SCOPES, execution_scope=SCOPE,
            terminal_reviewed=True, composed_full_review_complete=True, native_outputs_validated=True,
            primary_resources_replayed=True, shared_host_resources_reviewed=True,
            current_original_inventory_equality=False, original_ordinary_full_review_success=False,
            next_identity_authorized=False, downstream_adoption_complete=False, accuracy_evaluated=False,
            scientific_timings_admitted=False, uncontended_timing=False, automatic_retry=False,
            publication_ready=False, timing_disclosure=DISCLOSURE,
            whole_run_maximum_foreign_average_cores=environment["processes"]["maximum_observed_foreign_average_cores"],
            historical_failures_retained=[23986, 24033], observed_unix_ns=time.time_ns(),
            limitations=["Truthful new composed review, not success of an original failed reviewer.",
                "Semantic validation reuses independently successful unchanged original producer24031.",
                "Compatible conversion/scoring/history adoption is still required; no next identity authorized.",
                "Runtime brackets and periodic monitoring do not establish continuous integrity or isolation."]))
    except Exception as error:
        save(DESTINATION / "failure.json", dict(schema=SCHEMA, status="composed_review_failed",
            source=source, review_job_id=job, native_job_id=23985,
            error_type=type(error).__name__, error=str(error), terminal_reviewed=False,
            next_identity_authorized=False, automatic_retry=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(review(args.source_sha256), sort_keys=True))
