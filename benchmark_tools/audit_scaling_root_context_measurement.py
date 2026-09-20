"""Replay a frozen scaling task's measurement; do not admit timing or authorize runs."""

import argparse
import json
from pathlib import Path

from benchmark_tools.measure_root_context_scaling import classify_measurement
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.replay_scaling_root_context import replay
from benchmark_tools.verify_scaling_task_records import bind


def audit(plan_path, index, directory, recipe_path, recipe_sha, job):
    bound = bind(plan_path, index, directory, recipe_path, recipe_sha, job)
    replayed = replay(Path(directory) / "measurement", job, bound["measured_argv"])
    disposition = classify_measurement(bound["verification"], replayed, job)
    if disposition["corroborated_native_outcome"] != bound["native_outcome"]:
        raise ValueError("Bound task and replay do not corroborate the same native outcome")
    # Replay can be long: recheck records loaded before it, not only raw files.
    for item in [*bound["sources"], *bound["evidence"], bound["source"],
                 *replayed["evidence"], replayed["source"]]:
        check(item)
    return dict(status="scaling_task_measurement_audited", index=index, job_id=job,
        bound_task=bound, replay=replayed, disposition=disposition, source=record(__file__),
        scientific_timings_admitted=False, environmental_validity_established=False,
        native_outputs_validated=False, next_submission_authorized=False, publication_ready=False,
        limitations=["Composes task-record binding, raw replay and native-outcome classification only.",
            "Requires independent recipe archive, authorization, terminal scheduler/session and environmental audits.",
            "Native outputs or failure evidence still require validation; no automatic retry or timing admission."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "directory", "recipe", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--recipe-sha", required=True)
    parser.add_argument("--job", type=int, required=True)
    args = parser.parse_args()
    result = audit(args.plan, args.index, args.directory, args.recipe, args.recipe_sha, args.job)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
