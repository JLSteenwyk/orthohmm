import json
from pathlib import Path

import pytest

from benchmark_tools.probe_threadripper_native_fixture import fixture_run


@pytest.mark.parametrize("method", ["orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full"])
def test_fixture_relocation_preserves_scientific_arguments(tmp_path, method):
    plan = json.loads((Path(__file__).resolve().parents[2] /
        "benchmark_tools/results/threadripper_scaling_commands_20260928.json").read_text())
    template = next(r for r in plan["runs"] if r["native_method"] == method)
    before = json.dumps(template, sort_keys=True)
    source = tmp_path / "input"
    source.mkdir()
    for i in range(4):
        (source / f"S{i}.fa").write_text("".join(f">S{i}_{j}\nACDEFGHIK\n" for j in range(4)))
    run = fixture_run(template, source, tmp_path / "fixture_00",
                      Path("/dev/shm/orthohmm_fixture_unit_unused"))
    assert json.dumps(template, sort_keys=True) == before
    assert run["dataset"]["proteins"] == 16
    assert run["diagnostic_only"] is True
    assert not Path(run["prepared_input_directory"]).exists()
    for original, changed in zip(template["native_argv"], run["native_argv"]):
        if original != changed:
            assert original == template["prepared_input_directory"] or original.startswith(
                str(Path(template["measurement_directory"]).parent) + "/")
    (source / "S0.fa").write_text(">only\nACDE\n")
    with pytest.raises(ValueError, match="16-gene"):
        fixture_run(template, source, tmp_path / "fixture_00", Path("/dev/shm/unused"))
