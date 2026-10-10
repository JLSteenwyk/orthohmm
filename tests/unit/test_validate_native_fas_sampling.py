from collections import Counter
from decimal import localcontext
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import validate_native_fas_sampling as current


@pytest.fixture(autouse=True)
def caches_and_precision():
    current.checked_count_bounds.cache_clear()
    with localcontext() as ctx:
        ctx.prec=80
        yield
    current.checked_count_bounds.cache_clear()


def test_category_multiplicities_equal_actual_equiprobable_subsets():
    categories=[(2,None),(2,".25"),(1,".75")]
    genes=[(i,score) for i,(n,score) in enumerate(categories) for _ in range(n)]
    for c in range(6):
        counts=Counter(tuple(sum(i==index for i,_ in subset) for index in range(3))
                       for subset in combinations(genes,c))
        outcomes=current.finite_outcomes(categories,c)
        assert {s["selected"]:s["probability"] for s in outcomes}=={
            key:Fraction(value,sum(counts.values())) for key,value in counts.items()}
        assert sum(s["probability"] for s in outcomes)==1


def test_frozen_plan_contains_partial_and_literal_native_caps_and_unequal_strata():
    plan=current.plan()
    assert len(plan)==15 and len({x["name"] for x in plan})==15
    assert [x["name"] for x in plan if x["c"]==9000]==[
        "native_cap_sparse","native_cap_dense","native_cap_no_returns","native_cap_all_returns"]
    for cell in plan:
        P=sum(n for n,_ in cell["pre"]);M=sum(n for n,_ in cell["missing"])
        assert cell["k"]==min(P,max(1,round(cell["c"]*P/M)))
    assert current.ERROR==Fraction(1,320) and current.ALPHA==Fraction(1,20)


def test_independent_ratio_average_rejects_expected_count_substitution():
    row,states,target=current.validate_cell(current.plan()[4])
    assert target==Fraction(107,240)
    assert Fraction(row["native_ratio_average_fraction"])==target
    assert Fraction(row["expected_count_plugin_fraction"])==Fraction(9,20)
    assert Fraction(row["target_coverage_fraction"])>=1-2*current.ERROR
    assert sum(s["probability"] for s in states)==1


@pytest.mark.parametrize("index",[0,1,2,3,8,14])
def test_boundary_and_score_dependent_cells(index):
    row,_,_=current.validate_cell(current.plan()[index])
    assert Fraction(row["count_coverage_fraction"])>=1-current.ERROR
    assert Fraction(row["mean_coverage_fraction"])>=1-current.ERROR
    assert row["numerical_max_error"]<=current.TOLERANCE


def test_exact_integer_endpoint_checks_refuse_wrong_count_interval(monkeypatch):
    monkeypatch.setattr(current.kernel,"count_interval",lambda *args:(0,0))
    with pytest.raises(ValueError,match="Rejected count endpoints"):
        current.checked_count_bounds(4,2,1)


def test_numerical_disagreement_not_hidden_by_reference_coverage(monkeypatch):
    monkeypatch.setattr(current.kernel,"expected_native_mean",lambda *args:0.)
    with pytest.raises(ValueError,match="Target coverage/numerics failed"):
        current.validate_cell(current.plan()[4])


def test_coupled_eight_method_all_difference_projection(monkeypatch):
    row,states,target=current.validate_cell(current.plan()[4])
    monkeypatch.setattr(current,"JOINT_INDICES",tuple(range(8)))
    result=current.validate_joint([row]*8,[states]*8,[target]*8)
    assert result["method_count"]==8 and result["differences"]==28
    assert Fraction(result["coupled_mean_coverage_fraction"])>=1-current.ALPHA
    assert Fraction(result["coupled_difference_coverage_fraction"])>=1-current.ALPHA


def test_context_dependent_returns_are_rejected_not_claimed_covered():
    control=current.context_control()
    assert control["status"]=="fixed_return_assumption_rejected"
    assert Fraction(control["actual_expected_ratio_fraction"])==Fraction(22,45)
    assert Fraction(control["singleton_fixed_return_target_fraction"])==Fraction(8,15)
    assert control["coverage_claimed"] is False and control["native_admission"] is False


def test_occupied_output_refused_before_validation(tmp_path,monkeypatch):
    monkeypatch.setattr(sys,"argv",["validate","--protocol","missing","--output",str(tmp_path)])
    monkeypatch.setattr(current,"run",lambda *args:pytest.fail("Validation must not run"))
    with pytest.raises(FileExistsError):
        current.main()


def test_frozen_source_protocol_bindings_and_changed_kernel_refusal(monkeypatch):
    root=Path(__file__).resolve().parents[2]
    protocol=root/"benchmark_tools/results/NATIVE_FAS_SAMPLING_DESIGN_PROTOCOL_20261010.md"
    assert current.fingerprint(current.kernel.__file__)["sha256"]==current.KERNEL_SHA
    assert current.fingerprint(protocol)["sha256"]==current.PROTOCOL_SHA
    monkeypatch.setattr(current,"KERNEL_SHA","changed")
    with pytest.raises(ValueError,match="Changed frozen kernel/protocol"):
        current.run(protocol)


def test_reports_contain_finite_fraction_probabilities():
    row,_,_=current.validate_cell(current.plan()[4])
    assert json.loads(json.dumps(row,allow_nan=False))==row
