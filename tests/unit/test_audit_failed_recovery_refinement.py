import hashlib
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools.audit_failed_recovery_refinement import coverage, record


@pytest.mark.parametrize("text,names,count,valid", [
    ("a b\nc\n", ["a", "b", "c"], 2, True),
    ("a b\n\nc\n", ["a", "b", "c"], 2, True),
    ("a a\nc\n", ["a", "b", "c"], 2, False),
    ("a b\nb c\n", ["a", "b", "c"], 2, False),
    ("a b\nd\n", ["a", "b", "c"], 2, False),
    ("a b\n", ["a", "b", "c"], 1, False),
    ("a b\nc\n", ["a", "b", "c"], 3, False),
    ("a b\n", ["a", "b", "b"], 1, False),
    ("", [], 0, False),
])
def test_standalone_partition_readback(tmp_path, text, names, count, valid):
    path = tmp_path / "partition.txt"
    path.write_text(text)
    if valid:
        assert coverage(path, names, count) == dict(genes=len(names), groups=count)
    else:
        with pytest.raises(ValueError):
            coverage(path, names, count)
    assert record(path) == dict(path=str(path), bytes=len(text), sha256=hashlib.sha256(text.encode()).hexdigest())


def test_cli_requires_no_site_packages():
    from benchmark_tools import audit_failed_recovery_refinement as module
    done = subprocess.run([sys.executable, "-I", "-S", str(Path(module.__file__).resolve()), "--help"],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
