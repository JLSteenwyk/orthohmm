import hashlib
import json

import pytest

from benchmark_tools.run_ob_dependency_replay import load_plan, prepare


def test_execution_plan_requires_exact_hash(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text(json.dumps({"arms": ["a", "b"]}))
    checksum = hashlib.sha256(path.read_bytes()).hexdigest()
    assert load_plan(path,checksum) == {"arms": ["a", "b"]}
    with pytest.raises(ValueError):
        load_plan(path,None)
    with pytest.raises(ValueError):
        load_plan(path,"0" * 64)
    path.write_text("{}")
    with pytest.raises(ValueError):
        load_plan(path,checksum)


def test_prepare_never_overwrites_existing_directory(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path,tmp_path)
