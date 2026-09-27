import copy
import json

import pytest

from benchmark_tools import admit_integrated_full_ob as module


def stages(tmp_path):
    expected, outcomes = [], []
    for i in range(8):
        name, command, environment = f"stage{i}", ["python", str(i)], dict(LANG="C.UTF-8")
        expected.append((name, command, environment))
        (tmp_path / (name + ".log")).write_text("completed")
        module.save(tmp_path / (name + "_started.json"), dict(command=command, environment=environment))
        outcome = dict(command=command, returncode=0, log=module.record(tmp_path / (name + ".log")))
        module.save(tmp_path / (name + "_finished.json"), outcome)
        outcomes.append(outcome)
    return expected, dict(outcomes=outcomes)


def test_all_eight_stages_bound(tmp_path):
    expected, complete = stages(tmp_path)
    assert len(module.verify_stages(tmp_path, complete, expected)) == 24


def test_missing_stage_rejected(tmp_path):
    expected, complete = stages(tmp_path)
    complete["outcomes"].pop()
    with pytest.raises(ValueError, match="all eight"):
        module.verify_stages(tmp_path, complete, expected)


def test_reordered_stage_rejected(tmp_path):
    expected, complete = stages(tmp_path)
    complete["outcomes"][0], complete["outcomes"][1] = complete["outcomes"][1], complete["outcomes"][0]
    with pytest.raises(ValueError, match="completion differs"):
        module.verify_stages(tmp_path, complete, expected)


def test_failed_stage_cannot_be_hidden(tmp_path):
    expected, complete = stages(tmp_path)
    (tmp_path / "stage0_failed.json").write_text("{}")
    with pytest.raises(ValueError, match="Failed stage"):
        module.verify_stages(tmp_path, complete, expected)


def test_changed_log_rejected(tmp_path):
    expected, complete = stages(tmp_path)
    (tmp_path / "stage2.log").write_text("different")
    with pytest.raises(ValueError, match="completion differs"):
        module.verify_stages(tmp_path, complete, expected)


def test_changed_environment_rejected(tmp_path):
    expected, complete = stages(tmp_path)
    path = tmp_path / "stage0_started.json"
    value = json.loads(path.read_text())
    value["environment"]["LANG"] = "different"
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="environment differs"):
        module.verify_stages(tmp_path, complete, expected)


def test_no_admission_of_live_job(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|NodeList|ReqMem\n22337|RUNNING|0:0|00:01:00|32|bizon|128G\n")
    with pytest.raises(ValueError, match="COMPLETED"):
        module.audit(tmp_path, 22337, tmp_path / "admission")
    assert not (tmp_path / "admission").exists()


@pytest.mark.parametrize("cpus,node,memory", [("16", "bizon", "128G"), ("32", "other", "128G"), ("32", "bizon", "128Gc")])
def test_wrong_resources_rejected(tmp_path, monkeypatch, cpus, node, memory):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|NodeList|ReqMem\n22337|COMPLETED|0:0|01:00:00|{cpus}|{node}|{memory}\n")
    with pytest.raises(ValueError, match="resource allocation"):
        module.admission(tmp_path, 22337)


def test_wrong_job_rejected(tmp_path):
    with pytest.raises(ValueError, match="prespecified job"):
        module.admission(tmp_path, 22326)


def private_data(tmp_path):
    original = dict(dataset="installation_fixture", genes=16, fasta=[], references=[], uncertain=[])
    private = copy.deepcopy(original)
    for role, count in (("fasta", 4), ("references", 3), ("uncertain", 0)):
        folder = tmp_path / ("input" if role == "fasta" else "scoring_inputs/" + role)
        folder.mkdir(parents=True)
        for i in range(count):
            source = tmp_path / f"{role}{i}.txt"
            source.write_text(str(i))
            target = folder / source.name
            target.write_bytes(source.read_bytes())
            original[role].append(module.record(source))
            private[role].append(module.record(target))
    return original, private


def test_private_copy_byte_pins_bound(tmp_path):
    original, private = private_data(tmp_path)
    module.verify_private_inputs(original, private, tmp_path)


def test_extra_private_reference_rejected(tmp_path):
    original, private = private_data(tmp_path)
    (tmp_path / "scoring_inputs/references/extra.txt").write_text("extra")
    with pytest.raises(ValueError, match="Private input inventory"):
        module.verify_private_inputs(original, private, tmp_path)


def test_private_metadata_change_rejected(tmp_path):
    original, private = private_data(tmp_path)
    private["unexpected"] = True
    with pytest.raises(ValueError, match="metadata differs"):
        module.verify_private_inputs(original, private, tmp_path)
