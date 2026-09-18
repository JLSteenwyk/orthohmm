import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools.sequence_graph_evidence import sequence_evidence, sequence_settings
from benchmark_tools.validate_checked_replay_payload import validate


@pytest.mark.parametrize("stage,index", [("initial", 0), ("multipass", 1)])
def test_two_profile_off_stages(stage, index):
    metadata = dict(cpm_resolution=.1, seed=4, include_isolates=True, output_directory="/result/replay")
    sequence_settings(dict(stage=stage, index=index, output_directory="/result/replay"), metadata,
                      dict(output_root="/result"))


@pytest.mark.parametrize("field,value", [("stage", "profile_base"), ("stage", "profile_expanded"),
                                        ("index", 1), ("output_directory", "/elsewhere")])
def test_wrong_stage_rejected(field, value):
    manifest = dict(stage="initial", index=0, output_directory="/result/replay")
    manifest[field] = value
    with pytest.raises(ValueError):
        sequence_settings(manifest, dict(cpm_resolution=.1, seed=4, include_isolates=True,
            output_directory="/result/replay"), dict(output_root="/result"))


def test_plan_hash_before_parsing(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text("not json")
    with pytest.raises(ValueError, match="hash"):
        sequence_evidence(path, "wrong", "all_hits")


@pytest.mark.parametrize("kwargs", [dict(sequence_plan={}), dict(sequence_variant="all_hits"),
    dict(corrected_plan={}, sequence_plan={}, sequence_variant="all_hits")])
def test_provenance_modes_exclusive(tmp_path, kwargs):
    with pytest.raises(ValueError, match="provenance mode"):
        validate(tmp_path, {}, tmp_path, tmp_path, **kwargs)


@pytest.mark.parametrize("args", [["--sequence-plan", "plan"],
    ["--sequence-plan", "plan", "--sequence-plan-sha256", "hash", "--sequence-variant", "all_hits",
     "--corrected-plan", "hmm", "--corrected-plan-sha256", "hash"]])
def test_worker_cli_rejects_incomplete_or_mixed_modes(args):
    worker = Path(__file__).resolve().parents[2] / "benchmark_tools/checked_replay_payload_worker.py"
    result = subprocess.run([sys.executable, str(worker), "--root", "/missing", "--payload", "/missing",
        "--manifest", "/missing", "--manifest-sha256", "hash", *args], capture_output=True, text=True)
    assert result.returncode == 2
    assert "error:" in result.stderr
