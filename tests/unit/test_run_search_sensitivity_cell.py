import json
import os
from pathlib import Path
import sys
import time

import pytest

from benchmark_tools import run_search_sensitivity_cell as worker


def test_diamond_commands_are_species_specific(tmp_path):
    build, search = worker.diamond_commands("/diamond", tmp_path / "all.fasta",
                                           tmp_path / "species.fasta", tmp_path)
    assert build[1] == "makedb" and build[build.index("--in") + 1].endswith("species.fasta")
    assert search[search.index("--query") + 1].endswith("all.fasta")
    assert search[search.index("--db") + 1].endswith("species.dmnd")
    for option, value in {"--evalue": "1.0", "--threads": "4", "--matrix": "BLOSUM62",
                          "--gapopen": "11", "--gapextend": "1", "--comp-based-stats": "1",
                          "--masking": "1", "--max-target-seqs": "0", "--max-hsps": "1"}.items():
        assert search[search.index(option) + 1] == value
    assert "--very-sensitive" in search
    assert search[search.index("--outfmt") + 1:search.index("--very-sensitive")] == [
        "6", "qseqid", "sseqid", "qlen", "slen", "score", "bitscore", "evalue"]


@pytest.mark.parametrize("code", [0, 7])
def test_execute_preserves_exit_logs(tmp_path, code):
    stages = []
    def execute():
        worker.execute([sys.executable, "-c", f"print('test'); raise SystemExit({code})"],
                       tmp_path, "stage", dict(os.environ), time.monotonic() + 10, stages)
    if code:
        with pytest.raises(RuntimeError, match="exited 7"):
            execute()
    else:
        execute()
    assert stages[0]["returncode"] == code
    assert "test" in (tmp_path / "stage.log").read_text()
    assert (tmp_path / "stage.time.txt").is_file()


def test_deadline_does_not_launch(tmp_path):
    stages = []
    with pytest.raises(TimeoutError):
        worker.execute(["/does/not/exist"], tmp_path, "stage", {}, time.monotonic() - 1, stages)
    assert len(stages) == 1 and not list(tmp_path.iterdir())


def test_timeout_records_killed_child(tmp_path):
    stages = []
    with pytest.raises(worker.subprocess.TimeoutExpired):
        worker.execute([sys.executable, "-c", "import time; time.sleep(10)"],
                       tmp_path, "stage", dict(os.environ), time.monotonic() + .1, stages)
    assert stages[0]["returncode"] == -9


@pytest.fixture
def plan(tmp_path, monkeypatch):
    source = tmp_path / "input"
    source.mkdir()
    (source / "a.fasta").write_text(">gene\nACDEFG\n")
    row = dict(input=str(source), inputs=[worker.record(source / "a.fasta")], species=1,
               genes=1, seed=1, condition="test", split="calibration")
    staging = tmp_path / "runtime/staging.json"
    staging.parent.mkdir()
    staging.write_text("{}")
    python = staging.parent / "venv_clean/bin/python"
    python.parent.mkdir(parents=True)
    python.symlink_to(sys.executable)
    install = tmp_path / "install.json"
    install.write_text(json.dumps(dict(checked_records=[worker.record(staging)])))
    monkeypatch.setattr(worker, "INSTALL_SHA", worker.record(install)["sha256"])
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(dict(datasets=[row] * 70, diamond=worker.record(sys.executable),
                                   checked_evidence=[worker.record(install)])))
    monkeypatch.setattr(worker, "PLAN_SHA", worker.record(path)["sha256"])
    return path


def test_success_and_no_retry(plan, tmp_path, monkeypatch):
    calls = []
    def fake_execute(command, directory, name, env, deadline, stages):
        calls.append(command)
        if name == "hmm":
            dest = Path(command[command.index("--output") + 1])
            dest.mkdir()
            (dest / "hits.tsv").write_text("")
            (dest / "receipt.json").write_text("{}")
            assert "-I" in command and "PYTHONPATH" not in env
        elif command[1] == "makedb":
            Path(command[command.index("--db") + 1]).write_bytes(b"db")
        else:
            Path(command[command.index("--out") + 1]).write_text("")
    monkeypatch.setattr(worker, "execute", fake_execute)
    output = tmp_path / "run"
    worker.run(plan, 0, output)
    receipt = json.loads((output / "execution.json").read_text())
    assert receipt["status"] == "native_completed_pending_independent_readback"
    assert len(calls) == 3 and len(receipt["outputs"]) == 5
    with pytest.raises(FileExistsError):
        worker.run(plan, 0, output)
    assert len(calls) == 3


def test_failure_receipt_stops_later_stages(plan, tmp_path, monkeypatch):
    def fail(*args):
        raise RuntimeError("inference failure")
    monkeypatch.setattr(worker, "execute", fail)
    with pytest.raises(RuntimeError, match="inference failure"):
        worker.run(plan, 0, tmp_path / "run")
    receipt = json.loads((tmp_path / "run/execution.json").read_text())
    assert receipt["status"] == "failed" and receipt["attempt"] == 1
    assert "inference failure" in receipt["error"]
    assert not (tmp_path / "run/diamond").exists()


def test_plan_tamper_and_index_fail_before_output(plan, tmp_path):
    with pytest.raises(ValueError, match="index"):
        worker.run(plan, -1, tmp_path / "run")
    plan.write_text(plan.read_text() + " ")
    with pytest.raises(ValueError, match="Unrecognized"):
        worker.run(plan, 0, tmp_path / "run")
    assert not (tmp_path / "run").exists()


def test_input_change_is_failure(plan, tmp_path):
    (tmp_path / "input/a.fasta").write_text(">changed\nACD\n")
    with pytest.raises(ValueError, match="inventory"):
        worker.run(plan, 0, tmp_path / "run")
    assert json.loads((tmp_path / "run/execution.json").read_text())["status"] == "failed"
