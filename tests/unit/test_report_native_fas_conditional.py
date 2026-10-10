from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import report_native_fas_conditional as current


def inputs():
    samples, populations = [], []
    for method in current.audit.METHODS:
        mean = (180*.6+8998*.4)/(180+8998)
        samples.append(dict(method=method,strata=dict(
            precomputed=dict(population=200,requested=180,saved=180,mean=.6),
            missing=dict(population=10000,requested=9000,saved=8998,mean=.4)),
            saved_mean=mean,omitted_requested_new_scores=2,intended_sample_size=9180))
        populations.append(dict(method=method,precomputed=200,missing=10000,eligible_pairs=10200,
                                native_counts_match=True,saved_lookup_strata_and_values_match=True,
                                precomputed_mean=.61,precomputed_score_sum=122.,database_historically_hash_bound=False))
    return (dict(status="fas_requested_sample_attrition_bounded",methods=samples,
                 benchmark_scores_changed=False,uncertainty_admitted=False),
            dict(status="retained_fas_eligible_populations_completed_with_reuse",methods=populations,
                 benchmark_scores_changed=False,uncertainty_admitted=False,publication_ready=False,
                 reuse=dict(historical_parser_hash_identity_established=False)))


def interval(*args):
    return dict(expected_native_mean_bounds=[.3,.5],success_count_bounds=[9990,10000],
                return_mean_bounds=[.38,.42],target="expected_native_post_attrition_ratio_under_fixed_design",
                component_error=.05/16)


def test_full_population_mean_and_native_sample_means_used_without_changing_observed_scores(monkeypatch):
    calls=[]
    def capture(*args):
        calls.append(args)
        return interval(*args)
    monkeypatch.setattr(current.kernel,"method_interval",capture)
    sample,population=inputs()
    original=deepcopy((sample,population))
    result=current.panel(sample,population)
    assert all(x==(10000,9000,180,.61,8998,.4,.05/16) for x in calls) and len(calls)==8
    assert (sample,population)==original
    assert [x["observed_native_mean"] for x in result["methods"]]==[x["saved_mean"] for x in sample["methods"]]
    assert result["all_eight_ranges_computed"] and len(result["contrasts"])==28
    assert all(x["zero_included"] for x in result["contrasts"])
    assert all(result[k] is False for k in ("observed_native_scores_changed","historical_scores_rerun",
                                           "unconditional_historical_interval_admission","biological_generalization_intervals",
                                           "other_endpoint_uncertainty_admitted","publication_ready"))
    text=current.render(result)
    assert "not biological error bars" in text and "Unknown G and mu_G" in text
    assert len([x for x in text.splitlines() if x.startswith("| ")])==40


def test_one_numerical_failure_keeps_original_score_and_all_unavailable_contrasts(monkeypatch):
    calls=[]
    def fail_once(*args):
        calls.append(args)
        if len(calls)==1:
            raise ValueError("Hypergeometric mass not normalized")
        return interval(*args)
    monkeypatch.setattr(current.kernel,"method_interval",fail_once)
    sample,population=inputs()
    result=current.panel(sample,population)
    assert not result["all_eight_ranges_computed"] and len(calls)==8
    row=result["methods"][0]
    assert row["interval"] is None and row["observed_native_mean"]==sample["methods"][0]["saved_mean"]
    assert row["error"]["message"]=="Hypergeometric mass not normalized"
    assert sum(x["conditional_expected_difference_bounds"] is None for x in result["contrasts"])==7
    assert "partial numerical report" in current.render(result) and "Unavailable" in current.render(result)
    json.dumps(result,allow_nan=False)


@pytest.mark.parametrize("change",["method_order","population_mean","sum","observed_mean","count"])
def test_changed_identity_or_scalar_mapping_refused_before_intervals(monkeypatch,change):
    sample,population=inputs()
    if change=="method_order": sample["methods"].reverse()
    elif change=="population_mean": population["methods"][0]["precomputed_mean"]=float("nan")
    elif change=="sum": population["methods"][0]["precomputed_score_sum"]=100.
    elif change=="observed_mean": sample["methods"][0]["saved_mean"]=.8
    else: sample["methods"][0]["strata"]["precomputed"]["requested"]+=1
    monkeypatch.setattr(current.kernel,"method_interval",lambda *args:pytest.fail("No interval on invalid input"))
    with pytest.raises(ValueError):
        current.panel(sample,population)


def prerequisite_inputs():
    root=Path(__file__).resolve().parents[2]
    return [json.loads((root/current.PINS[name][0]).read_text()) for name in ("validation","context")]


@pytest.mark.parametrize("change",["late_population_mean","parser_scope"])
def test_all_scalar_mappings_validated_before_first_interval(monkeypatch,change):
    sample,population=inputs()
    if change=="late_population_mean": population["methods"][-1]["precomputed_mean"]=float("nan")
    else: population["reuse"]["historical_parser_hash_identity_established"]=None
    monkeypatch.setattr(current.kernel,"method_interval",lambda *args:pytest.fail("No partially checked application"))
    with pytest.raises(ValueError): current.panel(sample,population)


@pytest.mark.parametrize("change",["none","validation_status","coverage","control","context","admission"])
def test_prerequisite_scope_checks(change):
    validation,context=prerequisite_inputs()
    if change=="validation_status": validation["status"]="unfinished"
    elif change=="coverage": validation["cells"][0]["target_coverage_fraction"]="1/2"
    elif change=="control": validation["context_control"]["coverage_claimed"]=True
    elif change=="context": context["assessment"]["variations"]=["changed score"]
    elif change=="admission": validation["native_sampling_law_admitted"]=True
    if change=="none": current.prerequisites(validation,context)
    else:
        with pytest.raises(ValueError): current.prerequisites(validation,context)


def test_changed_prospective_protocol_refused_before_panel(monkeypatch):
    root=Path(__file__).resolve().parents[2]
    protocol=root/"benchmark_tools/results/NATIVE_FAS_CONDITIONAL_REPORTING_PROTOCOL_20261010.md"
    assert current.fingerprint(protocol)["sha256"]==current.PROTOCOL_SHA
    monkeypatch.setattr(current,"PROTOCOL_SHA","changed")
    monkeypatch.setattr(current,"panel",lambda *args:pytest.fail("No native application on changed pin"))
    with pytest.raises(ValueError,match="Changed prospective reporting protocol"):
        current.run(root,protocol)


def test_changed_scipy_runtime_refused_before_native_application(monkeypatch):
    root=Path(__file__).resolve().parents[2]
    protocol=root/"benchmark_tools/results/NATIVE_FAS_CONDITIONAL_REPORTING_PROTOCOL_20261010.md"
    monkeypatch.setattr(current.scipy,"__version__","unvalidated")
    monkeypatch.setattr(current,"panel",lambda *args:pytest.fail("No native interval under unvalidated runtime"))
    with pytest.raises(ValueError,match="Validated numerical runtime differs"):
        current.run(root,protocol)


def test_occupied_output_refused_before_reporting(tmp_path,monkeypatch):
    monkeypatch.setattr(sys,"argv",["report","--protocol","missing","--output",str(tmp_path)])
    monkeypatch.setattr(current,"run",lambda *args:pytest.fail("No duplicate producer"))
    with pytest.raises(FileExistsError): current.main()


def test_fatal_input_failure_is_preserved_in_fresh_namespace(tmp_path,monkeypatch):
    output=tmp_path/"new"
    monkeypatch.setattr(sys,"argv",["report","--protocol","missing","--output",str(output)])
    def fail(*args): raise ValueError("Changed input identity")
    monkeypatch.setattr(current,"run",fail)
    assert current.main()==1
    record=json.loads((output/"failure.json").read_text())
    assert record["message"]=="Changed input identity" and record["publication_ready"] is False
    assert not (output/"report.json").exists()
