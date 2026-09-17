import json
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.fixture(scope="module")
def binary(tmp_path_factory):
    compiler = shutil.which("gcc")
    if compiler is None:
        pytest.skip("GCC required for native load fixture")
    source = Path(__file__).resolve().parents[2] / "benchmark_tools/collector_load_fixture.c"
    output = tmp_path_factory.mktemp("collector_load") / "load"
    subprocess.run([compiler, "-O3", "-fopenmp", str(source), "-o", str(output)], check=True)
    return output


def test_deterministic_workers(binary):
    outputs = [subprocess.check_output([str(binary), "1000", "4"]) for _ in range(2)]
    assert outputs[0] == outputs[1]
    result = json.loads(outputs[0])
    assert result["workers"] == 4
    assert result["iterations_per_worker"] == 1000
    assert len(result["checksums"]) == 4
    one = json.loads(subprocess.check_output([str(binary), "1000", "1"]))
    assert one["checksums"] == result["checksums"][:1]


def test_dynamic_chunks_and_partial_tail(binary):
    a = subprocess.check_output([str(binary), "10000001", "2"])
    b = subprocess.check_output([str(binary), "10000001", "2"])
    assert a == b
    one = json.loads(subprocess.check_output([str(binary), "10000001", "1"]))
    assert one["checksums"] == json.loads(a)["checksums"][:1]


@pytest.mark.parametrize("args", [[], ["0", "4"], ["-1", "4"], ["x", "4"], ["100", "0"], ["100", "257"]])
def test_invalid_arguments(binary, args):
    assert subprocess.run([str(binary), *args]).returncode == 2
