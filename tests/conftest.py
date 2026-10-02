# global fixtures can go here
from pathlib import Path

import pytest


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


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")
    config.addinivalue_line("markers", "slow: mark as slow test")
