import pytest
from benchmark_tools.continue_full_qfo_cap400 import require_recovered_control


def scores():
    return {name: .8 for name in ("VGNC F", "SwissTrees F", "TreeFam-A F", "EC", "GO", "FAS")}


def test_control_gate_requires_all_six_verified_scores():
    report = {"status": "cached_native_control_and_all_six_scores_verified", "scores": scores()}
    require_recovered_control(report)
    for status in ("native_recovered_scoring_pending", "failed_native_recovery", "native_recovery_running"):
        with pytest.raises(ValueError, match="independently verified"):
            require_recovered_control({**report, "status": status})
    partial = scores()
    del partial["GO"]
    with pytest.raises(ValueError, match="endpoint"):
        require_recovered_control({**report, "scores": partial})
