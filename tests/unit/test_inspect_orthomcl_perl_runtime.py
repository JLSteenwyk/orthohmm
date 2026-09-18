import pytest

from benchmark_tools import inspect_orthomcl_perl_runtime as module


def fixture(tmp_path):
    source = tmp_path / "lib/module.pm"
    source.parent.mkdir()
    source.write_text("fixture")
    row = {**module.record(source), "kind": "file"}
    runtime = {"records": [row], "external_symlinks": module.EXTERNAL.copy()}
    observed = {"versions": module.VERSIONS.copy(), "search_path": [str(source.parent)],
                "loaded_modules": {"module.pm": str(source)}, "mapped_files": [str(source)]}
    cwd = tmp_path / "probe"
    cwd.mkdir()
    return runtime, observed, cwd


def test_loaded_records_bound(tmp_path):
    runtime, observed, cwd = fixture(tmp_path)
    assert len(module.validate(runtime, observed, cwd)) == 1


def test_guarded_do_script_requires_explicit_helper_binding(tmp_path):
    runtime, observed, cwd = fixture(tmp_path)
    driver = tmp_path / "driver.pl"
    driver.write_text("1;")
    observed["loaded_modules"][str(driver)] = str(driver)
    with pytest.raises(ValueError, match="outside snapshot"):
        module.validate(runtime, observed, cwd)
    helper = module.record(driver)
    assert len(module.validate(runtime, observed, cwd, [helper])) == 2
    driver.write_text("changed")
    with pytest.raises(ValueError):
        module.validate(runtime, observed, cwd, [helper])


@pytest.mark.parametrize("problem", ["version", "external", "relative", "hook", "missing", "changed", "empty", "invalid_path"])
def test_reject_runtime_drift(tmp_path, problem):
    runtime, observed, cwd = fixture(tmp_path)
    if problem == "version":
        observed["versions"]["perl"] = "other"
    elif problem == "external":
        runtime["external_symlinks"].append("/foreign")
    elif problem == "relative":
        observed["search_path"].append("../foreign")
    elif problem == "hook":
        observed["search_path"].append({"hook": True})
    elif problem == "missing":
        runtime["records"] = []
    elif problem == "changed":
        runtime["records"][0]["sha256"] = "changed"
    elif problem == "empty":
        observed["mapped_files"] = []
    else:
        observed["loaded_modules"]["other"] = {"hook": True}
    with pytest.raises(ValueError):
        module.validate(runtime, observed, cwd)


def test_dot_only_permitted_in_isolated_probe_directory(tmp_path):
    runtime, observed, cwd = fixture(tmp_path)
    observed["search_path"].append(".")
    (cwd / "probe.json").write_text("{}")
    (cwd / "probe.log").write_text("")
    module.validate(runtime, observed, cwd)
    (cwd / "foreign.pm").write_text("unreviewed")
    with pytest.raises(ValueError, match="current-directory"):
        module.validate(runtime, observed, cwd)


def test_probe_output_symlink_rejected(tmp_path):
    runtime, observed, cwd = fixture(tmp_path)
    observed["search_path"].append(".")
    (cwd / "probe.json").symlink_to(tmp_path / "lib/module.pm")
    with pytest.raises(ValueError, match="current-directory"):
        module.validate(runtime, observed, cwd)
