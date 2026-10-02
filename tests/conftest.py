# global fixtures can go here
from pathlib import Path
import shlex
import sys

import pytest

from tests.gnu_time_runtime import GnuTimeSubprocess, select_gnu_time


@pytest.fixture
def gnu_time_binary():
    try:
        return select_gnu_time()
    except FileNotFoundError as error:
        pytest.skip(str(error))


@pytest.fixture
def bind_test_gnu_time(gnu_time_binary, monkeypatch):
    """Bind only an explicitly requested smoke-test module's subprocess calls."""
    def bind(module):
        proxy = GnuTimeSubprocess(gnu_time_binary, module.subprocess)
        monkeypatch.setattr(module, "subprocess", proxy)
        return proxy
    return bind


@pytest.fixture
def synthetic_linux_boot_id(monkeypatch):
    """Stub only the boot counter for synthetic Linux workflow tests."""
    boot_path = Path("/proc/sys/kernel/random/boot_id")
    boot_id = "00000000-0000-4000-8000-000000000042"
    read_text = Path.read_text

    def synthetic_read(path, *args, **kwargs):
        if path == boot_path:
            return boot_id + "\n"
        return read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", synthetic_read)
    return boot_id


@pytest.fixture
def stage_historical_qfo_batch():
    """Relocate only ROOT and Python in a temporary handoff-test copy."""
    def stage(source, destination, root, *, python_calls=1):
        text = source.read_text()
        roots = [line for line in text.splitlines() if line.startswith("ROOT=")]
        historical_python = "/home/bizon/anaconda3/bin/python"
        if len(roots) != 1 or text.count(historical_python) != python_calls:
            raise ValueError("Historical batch relocation markers changed")
        relocated = text.replace(roots[0], "ROOT=" + shlex.quote(str(root)), 1)
        relocated = relocated.replace(historical_python, shlex.quote(sys.executable))
        with destination.open("x") as handle:
            handle.write(relocated)
        return destination
    return stage


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")
    config.addinivalue_line("markers", "slow: mark as slow test")
