import pytest

from benchmark_tools import run_three_kingdoms_matched_sonic as runner


@pytest.fixture
def plan(tmp_path, monkeypatch):
    directory = tmp_path / "inputs"
    directory.mkdir()
    inputs = []
    for index in range(12):
        path = directory / f"sp{index}.fasta"
        path.write_text(f">g{index}\nABC\n")
        inputs.append(runner.record(path))
    source = {"inputs": [{"files": {"staged": item}} for item in inputs],
              "counts": {"proteins": 443217}, "status": "retained_three_kingdoms_lineage_and_reference_verified"}
    monkeypatch.setattr(runner, "read_frozen", lambda *args: source)
    output = tmp_path / "new_output"
    return {"status": "matched_three_kingdoms_sonic_frozen_unrun", "output_root": str(output),
            "native_argv": runner.command(runner.ENTRYPOINT, output), "accuracy_admitted": False,
            "reuse_search": False, "resources": {"node": "bizon", "cpus": 32, "memory_gib": 192, "hours": 72},
            "input_source_audit": {"path": "fixture"}, "inputs": inputs}


def test_valid_frozen_plan(plan):
    assert str(runner.validate_plan(plan)) == plan["output_root"]


@pytest.mark.parametrize("mutation", ["mode", "threads", "resources", "reuse", "admission", "missing", "digest"])
def test_changed_plan_rejected(plan, mutation):
    if mutation == "mode":
        plan["native_argv"].extend(["-m", "fast"])
    elif mutation == "threads":
        plan["native_argv"][-1] = "64"
    elif mutation == "resources":
        plan["resources"]["memory_gib"] = 64
    elif mutation == "reuse":
        plan["reuse_search"] = True
    elif mutation == "admission":
        plan["accuracy_admitted"] = True
    elif mutation == "missing":
        plan["inputs"] = plan["inputs"][:-1]
    else:
        plan["inputs"] = [{**r, "sha256": "0" * 64} for r in plan["inputs"]]
    with pytest.raises(ValueError):
        runner.validate_plan(plan)


def test_unscheduled_execution_rejected_before_verification(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        runner.run(tmp_path / "absent", "0" * 64)


def test_check_only_creates_no_output(tmp_path, monkeypatch):
    output = tmp_path / "output"
    monkeypatch.setattr(runner, "verify", lambda *a: ({}, output, {}, {"test": "runtime"}))
    assert runner.run(tmp_path / "plan", "hash", True)["status"] == "preflight_passed_no_inference"
    assert not output.exists()


def test_existing_output_refused_even_for_preflight(tmp_path, monkeypatch):
    output = tmp_path / "output"
    output.mkdir()
    keep = output / "keep"
    keep.write_text("retained")
    monkeypatch.setattr(runner, "verify", lambda *a: ({}, output, {}, {}))
    with pytest.raises(FileExistsError):
        runner.run(tmp_path / "plan", "hash", True)
    assert keep.read_text() == "retained"
