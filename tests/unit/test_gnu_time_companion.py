import shutil
import subprocess
import sys

import pytest

from benchmark_tools.gnu_time_companion import command, parse


def test_parse_semantics():
    row = parse("elapsed_seconds\t1.20\nuser_seconds\t0.50\nsystem_seconds\t0.10\nmax_process_rss_kib\t1234\nexit_status\t0\n")
    assert row["elapsed_seconds"] == 1.2
    assert row["max_process_rss_kib"] == 1234
    assert "not simultaneous sum" in row["semantics"]["max_process_rss_kib"]


@pytest.mark.parametrize("bad", ["nan", "inf", "-1", "nope"])
def test_invalid_number(bad):
    with pytest.raises(ValueError):
        parse(f"elapsed_seconds\t{bad}\nuser_seconds\t0\nsystem_seconds\t0\nmax_process_rss_kib\t1\nexit_status\t0\n")


def test_reject_noise_and_missing_fields():
    with pytest.raises(ValueError):
        parse("Command exited with non-zero status 1\n")


@pytest.mark.parametrize("exit_code", [0, 7])
def test_actual_gnu_time_exit_and_output(tmp_path, exit_code):
    executable = shutil.which("time")
    if not executable:
        pytest.skip("External GNU time required")
    version = subprocess.check_output([executable, "--version"], text=True)
    if "GNU Time" not in version:
        pytest.skip("GNU time required")
    output = tmp_path / "accounting with spaces.tsv"
    argv = command([sys.executable, "-c", f"raise SystemExit({exit_code})"], output, executable)
    run = subprocess.run(argv, capture_output=True)
    assert run.returncode == exit_code
    assert parse(output.read_text())["exit_status"] == exit_code


def test_reject_relative_paths(tmp_path):
    with pytest.raises(ValueError):
        command(["echo", "x"], "relative.tsv")
    with pytest.raises(ValueError):
        command([], tmp_path / "out")
