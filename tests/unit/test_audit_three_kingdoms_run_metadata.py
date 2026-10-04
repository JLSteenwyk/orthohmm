import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import audit_three_kingdoms_run_metadata as audit
from benchmark_tools.assemble_benchmark_provenance import KEYS


def metrics(prediction):
    return dict(status="complete", command=["python", "__main__.py"], cwd="/recorded/cwd",
        wall_s=12, started_at_epoch_s=10, finished_at_epoch_s=22,
        user_cpu_s=24, system_cpu_s=3, peak_process_tree_rss_bytes=4096,
        rss_measurement="sampled tree RSS", harness=dict(exit_code=0,
            command=["python", "-m", "orthohmm"], git_commit="completion",
            git_dirty=True, output_manifest=[prediction], input_manifest=[], source_manifest=[]))


def test_metadata_events_preserve_recovery_history():
    result = audit.metadata_events("slurm_job_id\t1\nrecovery_slurm_job_id\t2\nrecovery_slurm_job_id\t3\n")
    assert result == [dict(key="slurm_job_id", value="1"),
                      dict(key="recovery_slurm_job_id", value="2"),
                      dict(key="recovery_slurm_job_id", value="3")]


@pytest.mark.parametrize("text", ["bad", "a\tb\tc", "a\t", "\tb", "\n"])
def test_invalid_metadata_events(text):
    with pytest.raises(ValueError):
        audit.metadata_events(text)


@pytest.mark.parametrize("damage", ["exit", "status", "output", "wall", "nan", "negative", "bool", "argv"])
def test_metrics_reject_unbound_or_invalid(damage):
    prediction = dict(path="/out", bytes=1, sha256="a" * 64)
    value = metrics(prediction)
    if damage == "exit":
        value["harness"]["exit_code"] = 1
    elif damage == "status":
        value["status"] = "running"
    elif damage == "output":
        value["harness"]["output_manifest"] = [dict(prediction, sha256="b" * 64)]
    elif damage == "wall":
        value["wall_s"] = 13
    elif damage == "nan":
        value["user_cpu_s"] = float("nan")
    elif damage == "negative":
        value["peak_process_tree_rss_bytes"] = -1
    elif damage == "bool":
        value["system_cpu_s"] = True
    else:
        value["command"] = "python __main__.py"
    with pytest.raises(ValueError):
        audit.metrics_details(value, prediction)


def test_metrics_retains_distinct_revisions_and_manifests():
    prediction = dict(path="/out", bytes=1, sha256="a" * 64)
    value = metrics(prediction)
    before = copy.deepcopy(value)
    result = audit.metrics_details(value, prediction)
    assert value == before
    assert result["harness_revision_at_completion"] == "completion"
    assert result["harness_git_dirty"]
    assert result["matched_output_manifest"] == [prediction]
    assert result["entrypoint_argv"] != result["harness_argv"]
    assert result["memory"]["unit"] == "bytes"


@pytest.fixture
def historical(tmp_path, monkeypatch):
    base = tmp_path / "benchmark_tools/results"
    base.mkdir(parents=True)
    rows, methods = [], []
    for key in KEYS:
        directory = tmp_path / key
        directory.mkdir()
        groups = directory / "orthogroups.txt"
        groups.write_text("OG0: gene1 gene2\n")
        pin = audit.record(groups)
        rows.append(dict(dataset="ThreeKingdoms", key=key, output_records=[pin],
                         declared_version="retained", resources=[]))
        methods.append(dict(key=key, result_dir=key,
                            provenance={"orthogroups_sha256":pin["sha256"]}))
        if key.startswith("orthohmm"):
            (directory / "metrics.json").write_text(json.dumps(metrics(pin)))
            (directory / "source_commit.txt").write_text("launch\n")
        elif key == "orthomcl_1_4":
            (directory / "time.log").write_text("")
            (directory / "run_metadata.tsv").write_text("recovery_slurm_job_id\t2\nrecovery_slurm_job_id\t3\n")
        else:
            name = "conversion_time.log" if key.endswith("sequence_only") else "time.log"
            (directory / name).write_text('Command being timed: "tool -c 32"\n'
                'User time (seconds): 24\nSystem time (seconds): 3\n'
                'Elapsed (wall clock) time (h:mm:ss or m:ss): 0:12.00\n'
                'Maximum resident set size (kbytes): 4096\nExit status: 0\n')
    # The replaced Sonic historical output must not be read or accidentally joined.
    methods[4]["provenance"]["orthogroups_sha256"] = "b" * 64
    (tmp_path / KEYS[4] / "time.log").write_text("invalid historical Sonic timing")
    reg = base / "all_benchmark_provenance_20261004_v2/register.json"
    reg.parent.mkdir()
    reg.write_text(json.dumps({"rows":rows}))
    parity = base / "three_kingdoms_parity_20260907.json"
    parity.write_text(json.dumps({"methods":methods}))
    monkeypatch.setattr(audit,"REGISTER_SHA",audit.record(reg)["sha256"])
    monkeypatch.setattr(audit,"PARITY_SHA",audit.record(parity)["sha256"])
    return tmp_path


def test_actual_shape_and_scopes(historical):
    output = historical / "supplement.json"
    result = audit.audit(historical, output)
    assert json.loads(output.read_text()) == result
    rows = {r["key"]:r for r in result["rows"]}
    assert len(rows) == 7 and KEYS[4] not in rows
    assert len(rows["orthomcl_1_4"]["metadata_events"]) == 2
    assert rows["orthomcl_1_4"]["empty_timing_logs"] == ["time.log"]
    high = rows[KEYS[0]]
    assert high["retained_text"]["source_commit.txt"] == "launch"
    assert high["metrics"]["harness_revision_at_completion"] == "completion"
    checkpoint = rows["orthofinder_3_1_5_sequence_only"]["timing_records"][0]
    assert checkpoint["scope"].startswith("checkpoint conversion only")
    assert not result["historical_consumption_proven"] and not result["publication_ready"]


@pytest.mark.parametrize("damage",["groups","source","failed_time"])
def test_changed_evidence_or_failed_time_refuses_before_output(historical, damage):
    output = historical / "supplement.json"
    if damage == "groups":
        (historical / KEYS[0] / "orthogroups.txt").write_text("changed")
    elif damage == "source":
        p=historical / "benchmark_tools/results/three_kingdoms_parity_20260907.json"
        p.write_text(p.read_text()+" ")
    else:
        p=historical / "proteinortho_6_3_6/time.log"
        p.write_text(p.read_text().replace("Exit status: 0","Exit status: 1"))
    with pytest.raises(ValueError):
        audit.audit(historical,output)
    assert not output.exists()


def test_preserve_existing_output(historical):
    output=historical / "supplement.json"
    output.write_text("existing")
    with pytest.raises(FileExistsError):
        audit.audit(historical,output)
    assert output.read_text()=="existing"
