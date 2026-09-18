import pytest

from benchmark_tools.validate_dgx_native_panel import validate_review


@pytest.mark.parametrize("problem", [None, "status", "failure", "missing", "index", "method", "repeat", "timing", "row_status"])
def test_complete_review_identity(problem):
    spec = {"runs": [{"index": i, "native_method": "orthohmm_high_sensitivity", "proteomes": 4, "repeat": i % 3}
                     for i in range(27)]}
    review = {"status": "retained_panel_review_complete_not_timing_admission", "review_failures": 0,
        "scientific_timings_admitted": 0, "controlled_workload_verified": False,
        "runs": [{"index": i, "method": "orthohmm_high_sensitivity", "proteomes": 4, "repeat": i % 3,
                  "status": "metadata_and_payload_hashes_verified_not_admitted"} for i in range(27)]}
    if problem == "status":
        review["status"] = "failed"
    elif problem == "failure":
        review["review_failures"] = 1
    elif problem == "missing":
        review["runs"].pop()
    elif problem == "index":
        review["runs"][2]["index"] = 1
    elif problem == "method":
        review["runs"][2]["method"] = "orthofinder_full"
    elif problem == "repeat":
        review["runs"][2]["repeat"] = 0
    elif problem == "timing":
        review["scientific_timings_admitted"] = 27
    elif problem == "row_status":
        review["runs"][2]["status"] = "failed"
    if problem:
        with pytest.raises(ValueError):
            validate_review(review, spec)
    else:
        validate_review(review, spec)
