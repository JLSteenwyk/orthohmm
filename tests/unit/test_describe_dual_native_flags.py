from copy import deepcopy

import pytest

from benchmark_tools import describe_dual_native_flags as module
from benchmark_tools.measure_native_dual_bracket_step import evaluate
from tests.unit.test_measure_native_dual_bracket_step import prepare
from tests.unit.test_measure_native_hierarchy_step import evidence


def test_phase_bounds_keep_all_uninstrumented_time_uncertain():
    bounds = module.phase_bounds(dict(stages={"search": {"wall_s": 10.}, "phylogeny": {"wall_s": 20.}}), 33.)
    assert bounds[1]["earliest_start_s"] < 10
    assert bounds[1]["latest_start_s"] > 13
    assert bounds[1]["earliest_end_s"] < 30
    assert module.certain_phase(bounds, 14, 29) == "phylogeny"
    assert module.certain_phase(bounds, 11, 12) is None
    assert module.certain_phase(bounds, 29, 31) is None


@pytest.mark.parametrize("stages,wall", [({}, 1), ({"unknown": {"wall_s": 1}}, 2),
    ({"search": {"wall_s": -1}}, 2), ({"search": {"wall_s": float("nan")}}, 2),
    ({"search": {"wall_s": 3}}, 2), ({"search": {"wall_s": 1}}, float("nan"))])
def test_invalid_phase_metrics_rejected(stages, wall):
    with pytest.raises(ValueError):
        module.phase_bounds(dict(stages=stages), wall)


@pytest.mark.parametrize("raw", ["", "processes -1", "processes 1 2", "processes 1\nprocesses 2", "processes nan"])
def test_invalid_proc_counter_rejected(raw):
    with pytest.raises(ValueError):
        module.counter(dict(raw=dict(proc_stat=raw)), "processes")


def test_valid_proc_counter():
    assert module.counter(dict(raw=dict(proc_stat="cpu 1 2\nprocesses 123\nctxt 456\n")), "processes") == 123


@pytest.fixture
def observed(evidence):
    points, done = prepare(evidence)
    for i, point in enumerate(points):
        for host in [*point["host"], point["hierarchy_host_after"]]:
            lines = [line for line in host["raw"]["proc_stat"].splitlines()
                     if line.split() and line.split()[0] not in {"processes", "ctxt"}]
            host["raw"]["proc_stat"] = "\n".join(lines) + f"\nprocesses {100+i*5}\nctxt {200+i*10}\n"
    screen = evaluate(points, done, 21816)
    measured = dict(points=points, screening=screen)
    run = dict(index=0, method="orthohmm_high_sensitivity", job_id=21816,
        started_ns=done["started_ns"], screening=deepcopy(screen),
        narrow_flagged_intervals=screen["narrow_flagged_intervals"][:])
    return run, measured


def test_complete_descriptive_join_does_not_assign_missing_phases(observed):
    result = module.describe(*observed, [])
    assert len(result["intervals"]) == 1
    row = result["intervals"][0]
    assert row["process_creations"] == 5
    assert row["context_switches"] == 10
    assert row["certain_phase"] is None
    assert result["all"]["certain_phases"] == {"unresolved": 1}
    assert result["by_certain_phase"]["unresolved"]["all"] == result["all"]
    assert result["flagged"]["count"] + result["unflagged"]["count"] == 1


@pytest.mark.parametrize("fault", ["flags", "screen", "counter", "inventory", "clock"])
def test_mismatched_descriptive_inputs_fail(observed, fault):
    run, measured = observed
    if fault == "flags":
        run["narrow_flagged_intervals"].append(999)
    elif fault == "screen":
        measured["screening"]["narrow_flagged_intervals"].append(999)
    elif fault == "counter":
        measured["points"][1]["hierarchy_host_after"]["raw"]["proc_stat"] = "processes 0\nctxt 0\n"
    elif fault == "inventory":
        measured["points"].pop()
    else:
        run["screening"]["narrow_intervals"][0]["wall_s"] += 1
        measured["screening"] = deepcopy(run["screening"])
    with pytest.raises(ValueError):
        module.describe(run, measured, [])
