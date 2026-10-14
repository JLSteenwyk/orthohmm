from copy import deepcopy
import pytest

from benchmark_tools.recover_qfo_native_phylogeny import require_native_crash, REVISION


def report():
    return {"status": "finished_with_explicit_gaps", "source_revision": REVISION,
            "arms": [{"label": "control100", "status": "failed", "inputs_sources_runtime_reverified": True}]}


def test_requires_known_native_signal_crash_not_success_or_scoring_failure():
    source = report()
    metrics = {"status": "failed", "harness": {"exit_code": -11}}
    require_native_crash(source, metrics)
    for exit_code in (0, 1, -9, 245):
        with pytest.raises(ValueError, match="signal-11"):
            require_native_crash(source, {"status": "failed", "harness": {"exit_code": exit_code}})


@pytest.mark.parametrize("field,value", [("status", "running"), ("source_revision", "wrong")])
def test_rejects_wrong_or_nonterminal_producer(field, value):
    source = report()
    source[field] = value
    with pytest.raises(ValueError): require_native_crash(source, {"status": "failed", "harness": {"exit_code": -11}})


def test_requires_preserved_postflight_source_input_runtime_proof():
    source = deepcopy(report())
    source["arms"][0]["inputs_sources_runtime_reverified"] = False
    with pytest.raises(ValueError): require_native_crash(source, {"status": "failed", "harness": {"exit_code": -11}})
