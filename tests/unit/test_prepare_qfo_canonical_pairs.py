import pytest

from benchmark_tools.prepare_qfo_canonical_pairs import EXPECTED_PAIRS, require_counts, verify_execution


def test_exact_counts():
    require_counts(EXPECTED_PAIRS, EXPECTED_PAIRS, EXPECTED_PAIRS)


@pytest.mark.parametrize("index,value", [(0,0), (1,5959534), (2,5959534), (0,5959535.0), (2,True)])
def test_count_or_mapping_change_rejected(index, value):
    counts = [EXPECTED_PAIRS] * 3
    counts[index] = value
    with pytest.raises(ValueError):
        require_counts(*counts)


@pytest.mark.parametrize("mutation", [None, "job_id", "status", "plan", "result"])
def test_execution_binding(mutation):
    plan, result = {"sha256": "plan"}, {"sha256": "result"}
    execution = dict(status="readback_complete", job_id="22334", plan=plan, result=result)
    if mutation:
        execution[mutation] = "wrong"
        with pytest.raises(ValueError):
            verify_execution(execution, plan, result)
    else:
        verify_execution(execution, plan, result)
