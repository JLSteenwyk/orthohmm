from copy import deepcopy
import json
import shutil

import pytest

from benchmark_tools import measure_native_root_context as module
from benchmark_tools.replay_native_root_context import replay
from tests.unit.test_measure_native_lineage_step import extend_lineage
from tests.unit.test_measure_native_hierarchy_step import evidence, setup_measure


def add_context(points):
    for index, point in enumerate(points):
        start = point["host"][1]["finished_monotonic_ns"] + 10
        ids = {name: [1, i+200] for i, name in enumerate(module.context.SCOPES)}

        def item(i, raw):
            return dict(started_ns=start+i*10, finished_ns=start+i*10+1, raw=raw)

        host = f"cpu {100+index} 0 100 500 0 0 0 0 0 0\n"
        point["root_context"] = dict(status="root_cpu_context_v1", ticks=point["ticks"],
            boot_before=point["lineage"]["boot_before"], boot_after=point["lineage"]["boot_before"],
            identities_before=ids, identities_after=deepcopy(ids), host_before=item(0, host),
            members_before=item(1, "2\n"), members_after=item(6, "2\n"), host_after=item(7, host),
            rows=[dict(item(i+2, f"usage_usec {100+index*(i+1)}\n"), scope=name)
                  for i, name in enumerate(module.context.SCOPES)])


def test_context_does_not_change_original_screens(evidence):
    points, done = extend_lineage(evidence)
    original = module.lineage.evaluate(points, done, 21816)
    add_context(points)
    result = module.evaluate(points, 21816)
    assert module.lineage.evaluate(points, done, 21816) == original
    assert result["intervals"][0]["root_minus_system_cpu_usec"] == -1
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["missing", "boot", "ticks", "early", "overlap", "counter"])
def test_invalid_supplementary_context_rejected(evidence, fault):
    points, done = extend_lineage(evidence)
    add_context(points)
    context = points[0]["root_context"]
    if fault == "missing":
        points[0].pop("root_context")
    elif fault == "boot":
        context["boot_before"] = context["boot_after"] = "other"
    elif fault == "ticks":
        context["ticks"] += 1
    elif fault == "early":
        context["host_before"]["started_ns"] = 0
    elif fault == "overlap":
        context["host_after"]["finished_ns"] = points[1]["host"][0]["started_monotonic_ns"]
    else:
        context["rows"][0]["raw"] = "usage_usec -1\n"
    with pytest.raises((ValueError, KeyError)):
        module.evaluate(points, 21816)


def test_wrapper_and_complete_replay(tmp_path, monkeypatch, evidence):
    import tests.unit.test_measure_native_hierarchy_step as fixture
    points, done = extend_lineage(evidence)
    add_context(points)
    monkeypatch.setattr(fixture, "module", module.lineage)
    directory = setup_measure(tmp_path, monkeypatch, evidence)
    memory = dict(scope=module.lineage.interval_point(points[-1], 21816)["native_cpu_scope"],
        errors=[], raw={"memory.current": "100", "memory.peak": "200"},
        started_ns=done["finished_ns"]+10, finished_ns=done["finished_ns"]+20)
    monkeypatch.setattr(module.lineage, "step_memory", lambda *args: memory)
    pending = iter(points)
    monkeypatch.setattr(module, "read_point", lambda *args: next(pending))
    measured, context = module.measure(["/usr/bin/true"], directory, 21816)
    result = replay(directory, 21816, ["/usr/bin/true"])
    assert result["context"] == context["context"]
    assert result["lineage"]["native_wall_s"] == measured["native_wall_s"]
    assert not result["scientific_timings_admitted"]
    relocated = tmp_path / "relocated"
    shutil.copytree(directory, relocated)
    assert replay(relocated, 21816, ["/usr/bin/true"])["context"] == result["context"]
    path = directory / "root_context_report.json"
    value = json.loads(path.read_text())
    value["context"]["intervals"][0]["root_minus_system_cpu_usec"] = 999
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="does not reproduce"):
        replay(directory, 21816, ["/usr/bin/true"])


def test_context_failure_retains_raw_partial(tmp_path, monkeypatch, evidence):
    points, _ = extend_lineage(evidence)
    monkeypatch.setattr(module.lineage, "read_point", lambda *args: points[0])

    def fail():
        raise module.context.RootContextError(ValueError("missing"), {"rows": []})

    monkeypatch.setattr(module.context, "snapshot", fail)
    with pytest.raises(module.context.RootContextError):
        module.read_point(123, "membership", 21816, tmp_path / "failed_point.json")
    value = json.loads((tmp_path / "failed_root_context.json").read_text())
    assert value["partial_context"]["rows"] == []
    assert not value["scientific_timings_admitted"]
