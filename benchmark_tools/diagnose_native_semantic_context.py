"""Probe the unchanged semantic kernel in its original caller context; no admission."""

import argparse
from collections import Counter
import json
from pathlib import Path
import traceback

from benchmark_tools.diagnose_native_partition import Bindings


def evaluate(context, validator):
    result = dict(status="semantic_probe_failed", context=context, semantic_result=None,
                  native_outputs_validated=False, accuracy_evaluated=False, terminal_reviewed=False,
                  next_identity_authorized=False, automatic_retry=False, resources_admitted=False,
                  historical_failure_cause_established=False)
    try:
        result["semantic_result"] = validator(context)
        result["status"] = "semantic_probe_passed"
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error),
                      traceback="".join(traceback.format_exception(type(error), error, error.__traceback__)))
        frames = []
        current = error.__traceback__
        while current is not None:
            frame = current.tb_frame
            owners = frame.f_locals.get("owners")
            item = dict(file=frame.f_code.co_filename, function=frame.f_code.co_name, line=current.tb_lineno)
            if isinstance(owners, dict):
                item["owner_count"] = len(owners)
                values = list(owners.values())
                if all(isinstance(value, str) for value in values):
                    item["owner_per_species"] = dict(sorted(Counter(values).items()))
                groups = frame.f_locals.get("value")
                if isinstance(groups, dict) and all(isinstance(g, list) for g in groups.values()):
                    genes = Counter(g for members in groups.values() for g in members)
                    item.update(parsed_unique_genes=len(genes), memberships=sum(genes.values()),
                                missing_count=len(owners.keys() - genes.keys()),
                                extra_count=len(genes.keys() - owners.keys()),
                                missing_sample=sorted(owners.keys() - genes.keys())[:20],
                                extra_sample=sorted(genes.keys() - owners.keys())[:20])
            frames.append(item)
            current = current.tb_next
        result["failure_frames"] = frames
    return result


def probe(diagnosis_ref):
    from benchmark_tools.native_factorial_allocated_execution import native_command
    from benchmark_tools.validate_native_factorial_outputs import validate_semantics
    import benchmark_tools.validate_native_factorial_outputs as kernel
    import benchmark_tools.validate_allocated_native_factorial_outputs as caller

    evidence = Bindings()
    diagnosis = evidence.read(diagnosis_ref["path"], diagnosis_ref)
    if (diagnosis.get("schema") != "native_partition_diagnosis_v1"
            or diagnosis.get("status") != "diagnosis_completed"
            or diagnosis.get("terminal_reviewed") is not False
            or diagnosis.get("native_outputs_validated") is not False):
        raise ValueError("Require the non-admitting partition diagnosis")
    for ref in diagnosis["evidence"]:
        evidence.bind(ref["path"], ref)
    failure = evidence.read(diagnosis["original_failure"]["path"], diagnosis["original_failure"])
    request_ref = failure["request"]
    request = evidence.read(request_ref["path"], request_ref)
    if request["index"] != diagnosis["index"] or request["job_id"] != diagnosis["job_id"]:
        raise ValueError("Original request identity differs")
    amendment_ref = request["amendment"]
    execution = evidence.read(amendment_ref["path"], amendment_ref)
    if request["plan"] != diagnosis["plan"] or execution["historical_plan"] != diagnosis["plan"]:
        raise ValueError("Original request/amendment plan differs")
    plan = evidence.read(diagnosis["plan"]["path"], diagnosis["plan"])
    baseline = evidence.read(plan["baseline"]["path"], plan["baseline"])
    run = plan["runs"][request["index"]]
    root = Path(run["output_root"])
    # Byte-for-byte semantic-context construction from the unchanged caller.
    context = dict(run, input_directory=str(root / "input"), cpu=32, threads_per_worker=4,
                   command=native_command(amendment_ref, run, baseline, metrics=True), cwd=baseline["core_root"],
                   aligner=baseline["tool_entrypoints"]["mafft"]["absolute_path"],
                   tree_builder=baseline["tool_entrypoints"]["FastTree"]["absolute_path"])
    for path in (__file__, kernel.__file__, caller.__file__):
        evidence.bind(Path(path).resolve())
    result = evaluate(context, validate_semantics)
    for ref in (result.get("semantic_result") or {}).get("checked_files", []):
        evidence.bind(ref["path"], ref)
    result.update(schema="native_semantic_context_diagnosis_v1", diagnosis=diagnosis_ref,
                  original_failure=diagnosis["original_failure"], request=request_ref,
                  index=request["index"], job_id=request["job_id"], evidence=evidence.finish(),
                  limitations=["The current kernel result does not determine the historical failure's cause.",
                               "No full allocated terminal review, resource replay, history authorization or scoring.",
                               "The nested kernel result retains its own non-admission scope."])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnosis", type=Path, required=True)
    parser.add_argument("--diagnosis-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink() or not args.output.parent.is_dir():
        raise ValueError("Require a fresh output in an existing directory")
    bindings = Bindings()
    bindings.bind(args.diagnosis)
    ref = bindings.files[str(args.diagnosis)]
    if ref["sha256"] != args.diagnosis_sha256:
        raise ValueError("Diagnosis checksum differs")
    result = probe(ref)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({k: result.get(k) for k in ("status", "index", "job_id", "error_type", "error",
                    "historical_failure_cause_established", "native_outputs_validated")}, sort_keys=True))


if __name__ == "__main__":
    main()
