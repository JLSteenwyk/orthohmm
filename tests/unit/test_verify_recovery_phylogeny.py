from pathlib import Path

import pytest

from benchmark_tools.verify_recovery_phylogeny import fixture_plan


def inputs():
    python = Path("/fresh/venv/bin/python")
    runtime = dict(executable=str(python), scientific_sources=[dict(
        path="/fresh/venv/lib/python3.10/site-packages/orthohmm/a.py")],
        dendropy_sources=[dict(path="/fresh/venv/lib/python3.10/site-packages/dendropy/a.py")])
    template = dict(fixture=True, cpu=1, attempts=1, checkpoint_reuse=False,
                    scoring=False, checked_records=[])
    return template, python, runtime


def test_fresh_scope_without_mutating_template():
    template, python, runtime = inputs()
    plan = fixture_plan(template, Path("/new"), python, runtime)
    assert plan["inference"] == "/new/inference"
    assert plan["python"] == str(python)
    assert len(plan["checked_records"]) == 2
    assert template["checked_records"] == []


@pytest.mark.parametrize("key,value", [("fixture", False), ("cpu", 32),
    ("attempts", 2), ("checkpoint_reuse", True), ("scoring", True)])
def test_reject_changed_scope(key, value):
    template, python, runtime = inputs()
    template[key] = value
    with pytest.raises(ValueError):
        fixture_plan(template, Path("/new"), python, runtime)


@pytest.mark.parametrize("fault", ["interpreter", "escape", "empty"])
def test_reject_wrong_installation(fault):
    template, python, runtime = inputs()
    if fault == "interpreter":
        runtime["executable"] = "/old/python"
    elif fault == "escape":
        runtime["scientific_sources"][0]["path"] = "/old/a.py"
    else:
        runtime["scientific_sources"] = []
    with pytest.raises(ValueError):
        fixture_plan(template, Path("/new"), python, runtime)
