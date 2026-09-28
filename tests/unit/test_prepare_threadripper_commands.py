import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_threadripper_commands import derive


def manifests():
    base = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    return [json.loads((base / name).read_text()) for name in
            ("publication_scaling_commands_20260917.json", "dgx_native_input_order_20260917.json")]


def scientific_tokens(argv, method):
    result = []
    i = 0
    while i < len(argv):
        if argv[i] in ("-f", "-o", "--metrics_json"):
            i += 2
        elif Path(argv[i]).is_absolute():
            i += 1
        else:
            result.append(argv[i])
            i += 1
    return result


def test_full_plan_preserves_science(tmp_path):
    plan, order = manifests()
    result = derive(plan, order, tmp_path / "output", Path("/dev/shm/test_threadripper_uncreated_root"))
    assert len(result["runs"]) == 27
    assert result["scientific_execution_authorized"] is False
    for old, new in zip(plan["runs"], result["runs"]):
        assert new["dataset"] == old["dataset"]
        assert new["cwd"] == old["cwd"]
        assert new["native_argv"][0] == old["native_argv"][0]
        assert scientific_tokens(new["native_argv"], new["native_method"]) == scientific_tokens(old["native_argv"], old["native_method"])
        assert new["native_argv"].count(new["prepared_input_directory"]) == 1
        assert Path(new["configuration"]["output"]).is_relative_to(tmp_path / "output")
        assert new["configuration"]["copy_inputs_to"] == new["prepared_input_directory"]
        if new["native_method"] == "orthofinder_full":
            assert new["native_argv"][-2:] == ["-o", new["configuration"]["output"]]
    assert len({r["prepared_input_directory"] for r in result["runs"]}) == 27
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("change", ["run_order", "bytes", "names", "disk_input", "tmpfs_output", "existing"])
def test_reject_protocol_drift(tmp_path, change):
    plan, order = manifests()
    out, inp = tmp_path / "output", Path("/dev/shm/test_threadripper_uncreated_root")
    if change == "run_order":
        plan["runs"].reverse()
    elif change == "bytes":
        order["datasets"][0]["inputs_in_native_order"][0]["sha256"] = "wrong"
    elif change == "names":
        order["datasets"][0]["native_order"].reverse()
    elif change == "disk_input":
        inp = tmp_path / "input"
    elif change == "tmpfs_output":
        out = Path("/dev/shm/test_threadripper_other_root")
    else:
        out.mkdir()
    with pytest.raises(ValueError):
        derive(plan, order, out, inp)
