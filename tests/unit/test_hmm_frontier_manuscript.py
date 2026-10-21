from copy import deepcopy
import math

import pytest

from benchmark_tools import build_hmm_frontier_manuscript as paper


def methods():
    return [dict(key="a", scores={"x": .7, "y": .9, "z": .2}),
            dict(key="b", scores={"x": .8, "y": .8, "z": .8}),
            dict(key="c", scores={"x": .6, "y": .8, "z": .9}),
            dict(key="tie", scores={"x": .7, "y": .9, "z": .2})]


def test_frontier_keeps_ties_and_records_all_dominators_without_mutation():
    data = methods(); original = deepcopy(data)
    rows = paper.frontier(data, ("x", "y"))
    assert data == original
    assert [r["key"] for r in rows if r["nondominated"]] == ["a", "b", "tie"]
    assert rows[2]["dominators"] == ["a", "b", "tie"]
    assert not paper.dominates((.7, .9), (.7, .9), (1, 1))


def test_pairwise_projection_does_not_determine_multivariate_membership():
    data = methods()
    assert not paper.frontier(data, ("x", "y"))[2]["nondominated"]
    assert paper.frontier(data, ("x", "y", "z"))[2]["nondominated"]


def test_unrounded_values_not_display_rounding_define_dominance():
    data = [dict(key="a", scores={"x": .7419645219079308}),
            dict(key="b", scores={"x": .7419468378090919})]
    assert round(data[0]["scores"]["x"], 3) == round(data[1]["scores"]["x"], 3)
    assert paper.frontier(data, ("x",))[1]["dominators"] == ["a"]


def test_directions_support_minimization_without_changing_accuracy_scopes():
    assert paper.dominates((.8, 10), (.7, 20), (1, -1))
    assert not paper.dominates((.8, 30), (.7, 20), (1, -1))
    assert paper.PRIMARY == ("OB", "VGNC", "SwissTrees")
    assert set(paper.ALL_F1) == {"OB", "VGNC", "SwissTrees", "TreeFam", "TK"}
    assert not set(paper.ALL_F1) & {"EC", "GO", "FAS"}


@pytest.mark.parametrize("value", [None, "NA", math.nan, math.inf, -math.inf])
def test_missing_or_nonfinite_objective_cannot_enter_frontier(value):
    data = methods(); data[0]["scores"]["x"] = value
    with pytest.raises(ValueError): paper.frontier(data, ("x", "y"))


@pytest.mark.parametrize("kind", ["duplicate_key", "duplicate_objective", "empty_objective", "missing_objective"])
def test_ambiguous_comparison_identity_is_rejected(kind):
    data = methods(); objectives = ("x", "y")
    if kind == "duplicate_key": data[0]["key"] = "b"
    elif kind == "duplicate_objective": objectives = ("x", "x")
    elif kind == "empty_objective": objectives = ()
    else: data[0]["scores"].pop("x")
    with pytest.raises(ValueError): paper.frontier(data, objectives)


@pytest.mark.parametrize("left,right,directions", [([], [], []), ([1], [1, 2], [1]),
    ([1], [1], [0]), ([1], [1], [2]), ([None], [1], [1])])
def test_invalid_dominance_dimensions_and_directions(left, right, directions):
    with pytest.raises(ValueError): paper.dominates(left, right, directions)


def test_uncertainty_preserves_saved_phase_intervals_without_resampling():
    values = dict(difference_percentage_points=.09, paired_percentile_ci=[-.05, .33],
                  bonferroni_percentile_ci=[-.08, .47])
    row = paper.uncertainty_row("weighted", values)
    assert row == dict(label="weighted", difference=.09, paired_percentile_ci=[-.05, .33],
                       bonferroni_percentile_ci=[-.08, .47])


@pytest.mark.parametrize("interval", [[1, -1], [math.nan, 1], [0], [0, 1, 2]])
def test_invalid_saved_uncertainty_is_not_drawn(interval):
    with pytest.raises(ValueError):
        paper.uncertainty_row("bad", dict(difference_percentage_points=0,
            paired_percentile_ci=interval, bonferroni_percentile_ci=[-1, 1]))


def test_existing_output_is_never_overwritten(tmp_path):
    marker = tmp_path / "original.txt"; marker.write_text("preserved")
    with pytest.raises(FileExistsError): paper.export({}, tmp_path)
    assert marker.read_text() == "preserved"


def test_partial_panel_is_not_exported_as_complete(tmp_path):
    output = tmp_path / "new"
    with pytest.raises(ValueError): paper.export(dict(methods=[]), output)
    assert not output.exists()


def test_changed_admission_fails_before_any_generation(tmp_path):
    p = tmp_path / paper.ARCHIVE; p.mkdir(parents=True)
    (p / "summary.json").write_text("{}")
    with pytest.raises(ValueError, match="Changed admitted input"):
        paper.collect(tmp_path)


def test_method_inventory_contains_all_fixed_configs_not_best_per_dataset():
    assert len(paper.METHODS) == 12
    keys = [m[0] for m in paper.METHODS]
    assert len(set(keys)) == 12
    assert "orthofinder_3_1_5_full" in keys and "orthofinder_3_1_5_sequence_only" in keys
    assert {"profile_weighting:none", "profile_weighting:henikoff", "profile_weighting:henikoff_cluster",
            "projection_recovery:projection_off"}.issubset(keys)


def test_manuscript_names_are_hmm2_without_renaming_historical_run_keys():
    own = paper.METHODS[:6]
    assert all(label.startswith("OrthoHMM2 ") for _, label, *_ in own)
    assert own[4][0] == "orthohmm_high_sensitivity"
    assert own[5][0] == "orthohmm_phylogeny_satellite_v2"


def timing_panel():
    runs, cells = [], []
    for method in ("sensitive", "phylogeny", "finder"):
        for size in (4, 8, 12):
            for repeat in range(3):
                index = len(runs); absent = index in (17, 20); eligible = index not in (0, 17, 20)
                runs.append(dict(index=index, method=method, proteomes=size, repeat=repeat,
                    comparative_timing_eligible=eligible,
                    status="reviewed_pre_native_abort" if absent else "reviewed_shared_observation" if eligible else "reviewed_excluded_attempt",
                    resources=None if absent else dict(wall_seconds=float(repeat + 1))))
            eligible_rows = [r for r in runs if r["method"] == method and r["proteomes"] == size and r["comparative_timing_eligible"]]
            cells.append(dict(method=method, proteomes=size, eligible_repeats=len(eligible_rows),
                resources=dict(wall_seconds=dict(median=2., minimum=1., maximum=3.) if len(eligible_rows) == 3 else None)))
    return dict(runs=runs, cells=cells)


def test_timing_keeps_aborts_missing_and_excluded_measurement_out_of_summaries():
    panel = timing_panel(); original = deepcopy(panel)
    paper.validate_timing(panel)
    assert panel == original
    assert panel["runs"][17]["resources"] is None
    assert panel["runs"][0]["resources"] is not None
    assert panel["cells"][0]["resources"]["wall_seconds"] is None


@pytest.mark.parametrize("kind", ["drop_abort", "zero_abort", "include_excluded", "partial_median", "wrong_median", "duplicate_identity"])
def test_invalid_timing_snapshots_are_not_presented_as_complete(kind):
    panel = timing_panel()
    if kind == "drop_abort": panel["runs"].pop(17)
    elif kind == "zero_abort": panel["runs"][17]["resources"] = dict(wall_seconds=0)
    elif kind == "include_excluded": panel["runs"][0]["comparative_timing_eligible"] = True
    elif kind == "partial_median": panel["cells"][0]["resources"]["wall_seconds"] = dict(median=2., minimum=1., maximum=3.)
    elif kind == "wrong_median": panel["cells"][1]["resources"]["wall_seconds"]["median"] = 1.
    else: panel["runs"][1]["index"] = 0
    with pytest.raises(ValueError): paper.validate_timing(panel)
