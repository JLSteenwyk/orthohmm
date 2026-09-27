import copy

import pytest

from benchmark_tools.admit_qfo_raw_tree_reuse import collect_records, verify_completion


def test_nested_records_deduplicated_and_normalized():
    r = dict(path="/x", bytes=3, sha256="abc")
    assert collect_records([r, dict(nested=[dict(r, label="extra")])]) == [r]


def test_conflicting_records_rejected():
    r = dict(path="/x", bytes=3, sha256="abc")
    with pytest.raises(ValueError):
        collect_records([r, dict(r, sha256="changed")])


def binding():
    plan, result_record, admission = [dict(path="/" + n, bytes=1, sha256=n) for n in ("plan", "result", "admission")]
    execution = dict(status="readback_complete", job_id="22330", plan=plan, result=result_record)
    result = dict(status="fresh_qfo_phylogeny_scientific_readback_complete", job_id=22329,
                  admission=admission, accuracy_evaluated=False, historical_scores_replaced=False,
                  partition=dict(label_invariant_equal=False))
    return execution, result, plan, result_record, admission


def test_valid_differences_do_not_block_reuse():
    verify_completion(*binding())


@pytest.mark.parametrize("index,key,value", [(0,"status","partial"), (0,"job_id","22329"),
    (0,"plan",{}), (0,"result",{}), (1,"status","partial"), (1,"job_id",22328),
    (1,"admission",{}), (1,"accuracy_evaluated",True), (1,"historical_scores_replaced",True)])
def test_completion_mutations_rejected(index, key, value):
    args = copy.deepcopy(binding())
    args[index][key] = value
    with pytest.raises(ValueError):
        verify_completion(*args)
