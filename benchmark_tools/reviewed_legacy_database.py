"""Bind fresh database evidence to the seven reviewed native O deletions."""

import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

REVIEW_SHA = "826e711979f0b262f31964e08b46bc39ff3cd621803074e1f71be213cc153a9e"


def validate_content(database, review):
    if (review["status"] != "reviewed_native_O_deletions_described_not_admitted"
            or database["status"] != "database_sequence_differences_require_review"
            or database["content"] != review["content"]):
        raise ValueError("Fresh database differs from reviewed native representation")
    content = database["content"]
    if (content["input_sequences"] != 984137 or content["database_sequences"] != 984137
            or content["exact_sequence_matches"] != 984130
            or content["input_residues"] != 440246934
            or content["database_residues"] != 440246927
            or content["exact_sequence_parity"] is not False
            or content["header_order_verified"] is not True
            or len(review["transformations"]) != 7
            or any(t["deleted_residue"] != "O" or len(t["positions_one_based"]) != 1
                   for t in review["transformations"])):
        raise ValueError("Unexpected reviewed representation scope")


def verify(root, database):
    path = root / "benchmark_tools/results/qfo_recovery_native_residue_deletions_20260925.json"
    identity = record(path)
    if identity["sha256"] != REVIEW_SHA:
        raise ValueError("Native representation review changed")
    review = json.loads(path.read_text())
    records = [identity, record(__file__), *review["checked_records"]]
    for item in records:
        check(item)
    validate_content(database, review)
    return dict(status="reviewed_native_representation_verified_not_exact_parity",
                exact_sequence_parity=False, transformations=review["transformations"],
                checked_records=records,
                limitations=["Seven native O deletions are retained, not corrected.",
                             "Query-side deletion is demonstrated by a native fixture, not a full query dump.",
                             "Reference exposure and counterfactual accuracy impact remain unevaluated."])
