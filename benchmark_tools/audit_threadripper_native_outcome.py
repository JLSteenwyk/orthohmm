"""Replay raw resource evidence and check successful native outputs independently."""

import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.replay_threadripper_scaling import replay
from benchmark_tools.validate_threadripper_outputs import validate


def audit(run, baseline, job_id, *, replay_fn=replay, validate_fn=validate):
    """Caller must bind run/baseline identity and separately audit runtime and host."""
    directory = Path(run["measurement_directory"])
    if not directory.is_absolute() or directory.resolve() != directory:
        raise ValueError("Require direct absolute measurement directory")
    paths = [directory / name for name in ("done.json", "command.json", "native.log")]
    if any(p.is_symlink() for p in paths):
        raise ValueError("Require direct native evidence files")
    evidence = [record(p) for p in paths]
    done, command = [json.loads(p.read_text()) for p in paths[:2]]
    if command.get("command") != run["native_argv"]:
        raise ValueError("Native command differs from expected run")
    reproduced = replay_fn(directory, job_id, run["native_argv"])
    if reproduced["measured"]["native"] != done:
        raise ValueError("Replayed native outcome differs from retained exit record")
    outcome = reproduced["native_outcome"]
    if outcome not in {"exited_zero", "exited_nonzero", "timed_out"}:
        raise ValueError("Unknown replayed native outcome")
    outputs = None
    if outcome == "exited_zero":
        outputs = validate_fn(run, dict(done, command=command["command"], cwd=run["cwd"]), baseline)
        if outputs.get("status") != "threadripper_native_outputs_checked":
            raise ValueError("Native output validation did not complete")
    for ref in evidence:
        check(ref)
    return dict(status="native_success_outputs_verified" if outputs else "native_failure_requires_review",
        job_id=job_id, native_outcome=outcome, replay=reproduced, outputs=outputs,
        evidence=evidence, source=record(__file__), automatic_retry=False,
        scientific_timings_admitted=False, next_submission_authorized=False,
        limitations=["Caller must bind frozen run/baseline/recipe, controller and runtime evidence separately.",
                     "Replayed failure is retained, not automatically classified as an eligible scientific failure.",
                     "Successful outputs and replay do not establish accuracy or quiet-host eligibility."])
