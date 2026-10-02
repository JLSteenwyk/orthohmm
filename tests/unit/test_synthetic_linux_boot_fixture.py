from pathlib import Path

import pytest


def test_boot_read_is_explicitly_synthetic(synthetic_linux_boot_id):
    assert synthetic_linux_boot_id == "00000000-0000-4000-8000-000000000042"
    assert Path("/proc/sys/kernel/random/boot_id").read_text().strip() == synthetic_linux_boot_id


def test_other_file_reads_are_not_stubbed(tmp_path, synthetic_linux_boot_id):
    path = tmp_path / "evidence.json"
    path.write_text('{"synthetic": true}\n')
    assert path.read_text(encoding="ascii") == '{"synthetic": true}\n'
    with pytest.raises(FileNotFoundError):
        (tmp_path / "missing.json").read_text()
