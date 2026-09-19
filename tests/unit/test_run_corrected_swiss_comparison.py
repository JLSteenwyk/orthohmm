import json

import pytest

from benchmark_tools import run_corrected_swiss_comparison as module
from benchmark_tools.audit_corrected_swiss_comparison import validate_inventory
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_bootstrap_corrected_swiss_comparators import fixture


@pytest.mark.parametrize("problem", [None, "order", "status", "count", "bool", "imputed", "conversion", "scores"])
def test_comparison_inventory_validation(problem):
    counts = fixture()
    methods = [{"key": row["method"], "status": "not_admitted", "scores": {"SwissTrees": None}}
               for row in counts["methods"]]
    comparison = {"status": "corrected_qfo_publication_comparison", "publication_ready": False,
                  "admitted_methods": 0, "methods": methods}
    if problem == "order":
        methods.reverse()
    elif problem == "status":
        methods[0]["status"] = "running"
    elif problem in ("count", "bool"):
        comparison["admitted_methods"] = 1 if problem == "count" else False
    elif problem in ("imputed", "conversion"):
        methods[0]["admission" if problem == "imputed" else "conversion"] = {}
    elif problem == "scores":
        methods[0]["scores"]["SwissTrees"] = 0
    if problem:
        with pytest.raises(ValueError):
            validate_inventory(comparison)
    else:
        validate_inventory(comparison)


def setup(tmp_path, monkeypatch):
    paths = [tmp_path / name for name in ("comparison.json", "baseline.json", "protocol.md", "release.md")]
    for path in paths:
        path.write_text("{}")
    monkeypatch.setattr(module, "PROTOCOL_SHA", record(paths[2])["sha256"])
    monkeypatch.setattr(module, "RELEASE_PROTOCOL_SHA", record(paths[3])["sha256"])
    counts = fixture()
    counts["checked_inputs"] = [record(path) for path in paths]
    monkeypatch.setattr(module, "audit", lambda *args: counts)
    return (paths[0], record(paths[0])["sha256"], paths[1], paths[2], paths[3], tmp_path / "result.json"), counts


@pytest.mark.parametrize("missing", [(), (6, 7), tuple(range(8))])
def test_real_kernel_integration_and_missing_status(tmp_path, monkeypatch, missing):
    args, counts = setup(tmp_path, monkeypatch)
    for index in missing:
        counts["methods"][index] = {"method": counts["methods"][index]["method"],
                                    "status": "not_admitted", "reason": "pending"}
    result = module.run(*args)
    assert result == json.loads(args[-1].read_text())
    assert result["estimated_contrasts"] == (8 if not missing else 6 if len(missing) == 2 else 0)
    assert result["complete_panel"] == (not missing)
    assert result["uncertainty_admitted"] == (len(missing) != 8)
    assert result["publication_ready"] is False


@pytest.mark.parametrize("problem", ["helper", "protocol", "release", "mutation", "controls", "audit"])
def test_fail_closed_before_output(tmp_path, monkeypatch, problem):
    args, counts = setup(tmp_path, monkeypatch)
    if problem == "helper":
        monkeypatch.setattr(module, "SOURCES", {"audit_corrected_swiss_comparison.py": "wrong"})
    elif problem in ("protocol", "release"):
        args[3 if problem == "protocol" else 4].write_text("changed")
    elif problem in ("mutation", "audit"):
        def audit(*unused):
            if problem == "audit":
                raise ValueError("Raw counts invalid")
            args[0].write_text("changed during analysis")
            return counts
        monkeypatch.setattr(module, "audit", audit)
    else:
        monkeypatch.setattr(module, "bootstrap", lambda *args: {"protocol_controls_match": False})
    with pytest.raises(ValueError):
        module.run(*args)
    assert not args[-1].exists()


@pytest.mark.parametrize("symlink", [False, True])
def test_existing_result_refused(tmp_path, monkeypatch, symlink):
    args, _ = setup(tmp_path, monkeypatch)
    if symlink:
        args[-1].symlink_to(tmp_path / "missing")
    else:
        args[-1].write_text("preserve")
    with pytest.raises(FileExistsError):
        module.run(*args)
