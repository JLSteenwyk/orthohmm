import hashlib
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import verify_orthobench_acquisition as module


def test_retained_inventory_complete():
    rows = module.expected_inputs(Path(__file__).resolve().parents[2] / "benchmark_tools/results")
    assert len(rows) == 93


@pytest.mark.parametrize("fault", [None, "changed", "missing", "extra", "symlink", "historical", "revision"])
def test_acquired_checkout_is_bound_to_git_and_input_bytes(tmp_path, monkeypatch, fault):
    paths = ["BENCHMARKS/Input/a.fa", "BENCHMARKS/RefOGs/RefOG001.txt",
             "BENCHMARKS/benchmark.py", "README.md"]
    expected = {}
    for name in paths:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"fixture\n")
        if name in paths[:2]:
            expected[name] = dict(bytes=8, sha256=hashlib.sha256(b"fixture\n").hexdigest())
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=tmp_path, stderr=subprocess.DEVNULL)
    git("init")
    git("add", ".")
    git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-m", "fixture")
    monkeypatch.setattr(module, "REVISION", git("rev-parse", "HEAD").decode().strip())
    monkeypatch.setattr(module, "expected_inputs", lambda results: expected)
    path = tmp_path / paths[0]
    if fault == "changed": path.write_text("modified")
    elif fault == "missing": path.unlink()
    elif fault == "extra": (path.parent / "extra.fa").write_text("extra")
    elif fault == "symlink":
        path.unlink()
        path.symlink_to(tmp_path / "README.md")
    elif fault == "historical": expected[paths[0]]["sha256"] = "0"*64
    elif fault == "revision": monkeypatch.setattr(module, "REVISION", "0"*40)
    if fault:
        with pytest.raises(ValueError): module.verify(tmp_path, tmp_path)
    else:
        result = module.verify(tmp_path, tmp_path)
        assert len(result["files"]) == 4
        assert result["redistribution_cleared"] is False
