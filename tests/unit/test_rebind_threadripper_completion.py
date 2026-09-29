import pytest
import json

from benchmark_tools import rebind_threadripper_completion as rebind
from benchmark_tools.rebind_threadripper_completion import changes
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_allowed_change():
    old = [dict(path="/a", sha256="old"), dict(path="/b", sha256="same")]
    new = [dict(path="/b", sha256="same"), dict(path="/a", sha256="new")]
    assert changes(old, new, {"/a"}) == [dict(path="/a", before=old[0], after=new[1])]


@pytest.mark.parametrize("new", [[], [dict(path="/b", sha256="old")],
    [dict(path="/a", sha256="old"), dict(path="/a", sha256="old")], [dict(path="/a", sha256="new")]])
def test_unexpected_inventory_or_content(new):
    with pytest.raises(ValueError):
        changes([dict(path="/a", sha256="old")], new, set())


def test_duplicate_prior():
    row = dict(path="/a", sha256="old")
    with pytest.raises(ValueError):
        changes([row, row], [row], {"/a"})


def test_failed_inventory_preserved(tmp_path, monkeypatch):
    prior_manifest = tmp_path / "manifest.json"
    prior_manifest.write_text(json.dumps(dict(roots=["/a"], records=[dict(path="/a", sha256="old")])))
    previous = tmp_path / "previous.json"
    previous.write_text(json.dumps(dict(baseline=record(prior_manifest), command_plan=record(prior_manifest),
        generator_source=record(prior_manifest), runtime_specs=[[str(prior_manifest), record(prior_manifest)["sha256"]]])))
    observed = dict(roots=["/a"], records=[dict(path="/b", sha256="new")])
    monkeypatch.setattr(rebind, "inventory", lambda _: observed)
    work, output = tmp_path / "work", tmp_path / "output.json"
    with pytest.raises(ValueError, match="inventory"):
        rebind.run(previous, record(previous)["sha256"], work, output)
    assert json.loads((work / "0.json").read_text()) == observed
    rejection = json.loads((work / "rejection.json").read_text())
    assert rejection["added"] == ["/b"] and rejection["removed"] == ["/a"]
    assert rejection["observed"] == record(work / "0.json")
    assert not output.exists()
