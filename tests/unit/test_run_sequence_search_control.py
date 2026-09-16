import os
import sys

import pytest

from benchmark_tools.run_sequence_search_control import run_phase


@pytest.mark.parametrize("exit_code", [0, 3])
def test_phase_records_process_outcome_and_resource_log(tmp_path, exit_code):
    command = [sys.executable, "-c", f"print('evidence'); raise SystemExit({exit_code})"]
    result = run_phase(command, tmp_path, "search", os.environ.copy())
    assert result["exit_code"] == exit_code
    assert result["argv"] == command
    assert result["wall_s"] >= 0
    assert "evidence" in (tmp_path / "search.log").read_text()
    assert f"Exit status: {exit_code}" in (tmp_path / "search.time.log").read_text()
    with pytest.raises(FileExistsError):
        run_phase(command, tmp_path, "search", os.environ.copy())
