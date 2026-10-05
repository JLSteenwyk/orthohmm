"""Run existing terminal review and frozen scoring once; never launch inference."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_native_factorial_attempt import review
from benchmark_tools.run_native_factorial_cost import ROOT, read, validate_plan, validate_request
from benchmark_tools.score_native_orthobench_attempt import score_attempt
from benchmark_tools.validate_native_factorial_outputs import require


def finish(request_ref, review_directory, score_directory):
    request = read(request_ref)
    validate_request(request, request["plan"], request["job_id"])
    run = validate_plan(read(request["plan"]))[request["index"]]
    require(run["dataset"] == "orthobench", "QfO requires its separate native endpoint workflow")
    destinations = [Path(p) for p in (review_directory, score_directory)]
    require(destinations[0] != destinations[1], "Review and score destinations must differ")
    for path in destinations:
        require(path.is_absolute() and path.resolve() == path and path.is_relative_to(ROOT),
                "Require separate direct repository destinations")
        require(not path.exists() and not path.is_symlink(), "Destination already exists")
    source = record(__file__)
    terminal_ref = review(request_ref, destinations[0])
    terminal = read(terminal_ref)
    require(terminal["status"] in {"native_success", "native_failure_retained"},
            "Unknown terminal review outcome")
    score_ref = None
    if terminal["status"] == "native_success":
        score_ref = score_attempt(request_ref, terminal_ref, destinations[1])
    check(source)
    result = {"schema": "native_orthobench_terminal_followup_v1", "request": request_ref,
              "source": source, "parent_native_job_id": request["job_id"],
              "review_score_scheduler_job_id": os.environ.get("SLURM_JOB_ID"),
              "terminal_review": terminal_ref, "score": score_ref,
              "status": "terminal_reviewed_and_scored" if score_ref else "native_failure_retained_unscored",
              "inference_reexecuted": False, "automatic_retry": False,
              "next_identity_released": False, "publication_ready": False}
    save(destinations[0] / "followup.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "review-directory", "score-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    require(sys.version_info[:2] == (3, 10), "Exact raw resource replay requires the matched Python3.10")
    require(record(__file__)["sha256"] == args.source_sha256, "Followup source checksum differs")
    ref = record(args.request)
    require(ref["sha256"] == args.request_sha256, "Followup request checksum differs")
    result = finish(ref, args.review_directory.absolute(), args.score_directory.absolute())
    print(json.dumps(result, sort_keys=True))
    if result["score"] is None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
