import pytest

from benchmark_tools.run_qfo_canonical_assessment import verify_stage


def fixture():
    identities = {k: dict(path="/" + k, bytes=1, sha256=k) for k in
                  ("pairs", "filtered_pairs", "mapping", "native_input", "plan", "result")}
    verified = dict(identities, status="canonical_conversion_independently_verified",
                    rows_compared=5959535, removed_mapping_pairs=0)
    stage = dict(identities, status="canonical_native_pairs_prepared_unscored", job_id="22335",
        participant="ohmm_qfo_canonical_20260927", semantics="native phylogenetically inferred pairs",
        accuracy_evaluated=False, total_pairs=5959535, retained_pairs=5959535, removed_mapping_pairs=0)
    execution = dict(status="conversion_complete", job_id="22335", plan=identities["plan"], result=identities["result"])
    return stage, verified, execution


def test_valid_binding():
    verify_stage(*fixture())


@pytest.mark.parametrize("index,key,value", [(0,"participant","other"), (0,"semantics","cliques"),
    (0,"accuracy_evaluated",True), (0,"job_id","22334"), (0,"filtered_pairs",{}),
    (0,"mapping",{}), (0,"native_input",{}), (0,"total_pairs",1),
    (1,"removed_mapping_pairs",1), (1,"status","partial"), (1,"rows_compared",1),
    (2,"job_id","wrong"), (2,"result",{}), (2,"status","partial")])
def test_invalid_binding(index, key, value):
    args = fixture()
    args[index][key] = value
    with pytest.raises(ValueError):
        verify_stage(*args)
