import copy

import pytest

from benchmark_tools.prepare_qfo_corrected_sonic import command, validate_runtime, DEFAULT_MODE, TOOLS


def runtime():
    return {"status": "read_only_current_resolution_observed", "version": "2.0.9",
            "default_mode": list(DEFAULT_MODE), "tools": {name: {"exit_code": 0} for name in TOOLS}}


def test_original_native_defaults_and_fresh_paths():
    assert command("/sonicparanoid", "/corrected/sonic") == [
        "/sonicparanoid", "-i", "/corrected/sonic/input", "-o", "/corrected/sonic/output", "-t", "32"]
    validate_runtime(runtime())


def test_relative_output_rejected():
    with pytest.raises(ValueError, match="absolute"):
        command("/sonicparanoid", "relative")


@pytest.mark.parametrize("key,value", [("version", "2.0.8"), ("status", "failed"),
                                      ("default_mode", ["fast", "diamond", "fast", 0.0])])
def test_changed_runtime_rejected(key, value):
    data = runtime()
    data[key] = value
    with pytest.raises(ValueError, match="Unexpected"):
        validate_runtime(data)


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_missing_or_failed_dependency(tool):
    data = runtime()
    missing = copy.deepcopy(data)
    del missing["tools"][tool]
    with pytest.raises(ValueError):
        validate_runtime(missing)
    data["tools"][tool]["exit_code"] = 1
    with pytest.raises(ValueError, match="probe failed"):
        validate_runtime(data)
