"""Replay both original lineage measurements and separately timed root context."""

import json
from pathlib import Path

from benchmark_tools.measure_native_root_context import evaluate, lineage_identity
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.replay_lineage_native_measurement import replay as replay_lineage


def replay(directory, job_id, command, expected_timeout_s=60):
    directory = Path(directory).absolute()
    report_record = record(directory / "root_context_report.json")
    retained = json.loads((directory / "root_context_report.json").read_text())
    original = replay_lineage(directory, job_id, command, expected_timeout_s)
    result = dict(status="native_root_context_measured", job_id=job_id,
        native_wall_s=original["native_wall_s"], context=evaluate(original["measured"]["points"], job_id),
        lineage_report=lineage_identity(directory),
        scientific_timings_admitted=False, environmental_validity_established=False)
    if retained != result:
        raise ValueError("Supplementary context report does not reproduce")
    check(report_record)
    for item in original["evidence"]:
        check(item)
    return dict(status="native_root_context_replayed", lineage=original, context=result["context"],
        evidence=[*original["evidence"], report_record], source=record(__file__),
        scientific_timings_admitted=False, environmental_validity_established=False)
