from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.prepare_qfo_sequence_graph import validate_payload, graph_command


def payload():
    variant = dict(status="frozen_rbnh_named_array_payload_bound", estimate=dict(
        genes=984137, species_slot_extent=78, graph_feasibility_admitted=False,
        total_peak_ram_bound_available=False, named_array_payload_lower_bytes=1024**3))
    return dict(status="admitted_qfo_graph_payload_estimated", source={"sha256": "source"},
                job_id="21798", variants={"all_hits": deepcopy(variant), "top100": deepcopy(variant)},
                graph_launched=False, graph_feasibility_admitted=False,
                accuracy_evaluated=False, publication_ready=False)


def validate(value, memory=None):
    validate_payload(value, {"sha256": "source"}, {"JobIDRaw": "21798"},
                     {"all_hits": 64, "top100": 32} if memory is None else memory)


def test_valid_review():
    validate(payload())


@pytest.mark.parametrize("memory", [{}, {"all_hits": 64}, {"all_hits": 1, "top100": 32},
                                     {"all_hits": True, "top100": 32}, {"all_hits": -1, "top100": 32}])
def test_invalid_review(memory):
    with pytest.raises(ValueError):
        validate(payload(), memory)


@pytest.mark.parametrize("field,value", [("genes", 251378), ("species_slot_extent", 12),
    ("graph_feasibility_admitted", True), ("total_peak_ram_bound_available", True)])
def test_wrong_estimate(field, value):
    report = payload()
    report["variants"]["all_hits"]["estimate"][field] = value
    with pytest.raises(ValueError):
        validate(report)


@pytest.mark.parametrize("field,value", [("job_id", "1"), ("source", {}), ("status", "failed"),
                                       ("accuracy_evaluated", True), ("graph_launched", True)])
def test_wrong_report(field, value):
    report = payload()
    report[field] = value
    with pytest.raises(ValueError):
        validate(report)


def test_native_command_preserves_profile_off_settings():
    command = graph_command(Path("/launcher"), Path("/checkpoint"), "hash", Path("/output"))
    assert "--fasta-directory" not in command
    for flag, value in (("--cpu", "32"), ("--matrix", "BLOSUM62"),
                        ("--cpm-resolution", "0.1"), ("--leiden-seed", "4"),
                        ("--checkpoint-sha256", "hash")):
        assert command[command.index(flag) + 1] == value
