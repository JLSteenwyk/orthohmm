import json
import sys

import pytest

from benchmark_tools.run_wgd_application import METHODS, check_copies, execute, pinned, select
from benchmark_tools.snapshot_orthohmm_input_order import record


def test_pinned_rejects_mutation(tmp_path):
    path = tmp_path / "manifest.json"
    path.write_text('{"a": 1}')
    identity = record(path)
    assert pinned(identity) == {"a": 1}
    path.write_text('{"a": 2}')
    with pytest.raises(ValueError, match="Changed pinned"):
        pinned(identity)


def test_select_requires_authorization_and_full_inventory(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text(json.dumps({"runs": [{"method": m} for m in METHODS], "output_root": str(tmp_path / "out")}))
    spec = {"command_plan": record(path), "execution_authorized": True, "purpose": "wgd_application"}
    assert select(spec, 3)[2] == tmp_path / "out/sonicparanoid"
    with pytest.raises(ValueError, match="index"):
        select(spec, -1)
    spec["execution_authorized"] = False
    with pytest.raises(ValueError, match="authorized"):
        select(spec, 0)


def test_copy_check_rejects_changes_and_extra_fasta(tmp_path):
    path = tmp_path / "a.fasta"
    path.write_text(">a\nMA\n")
    original = record(path)
    assert len(check_copies(tmp_path, [original])) == 1
    (tmp_path / "b.fasta").write_text(">b\nMA\n")
    with pytest.raises(ValueError, match="copies"):
        check_copies(tmp_path, [original])


@pytest.mark.parametrize("code", [0, 3])
def test_execute_preserves_exit_and_log(tmp_path, code):
    log = tmp_path / "native.log"
    result = execute([sys.executable, "-c", f"print('native'); raise SystemExit({code})"], tmp_path, log, 5)
    assert result["exit_code"] == code and not result["timed_out"]
    assert log.read_text() == "native\n"
    with pytest.raises(FileExistsError):
        execute([sys.executable, "-c", "pass"], tmp_path, log, 5)


def test_timeout_is_failure(tmp_path):
    result = execute([sys.executable, "-c", "import time; time.sleep(30)"], tmp_path, tmp_path / "native.log", 0.05)
    assert result["timed_out"] and result["exit_code"] != 0
