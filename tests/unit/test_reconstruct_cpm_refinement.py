import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import reconstruct_cpm_refinement as module
from benchmark_tools.audit_historical_profile_ablation import read_partition


@pytest.mark.parametrize("problem", [None, "changed", "duplicate", "missing", "unknown", "overwrite", "universe"])
def test_reconstruction_contract(tmp_path, problem):
    partition, expected, output = [tmp_path / name for name in ("input", "expected", "output")]
    partition.write_text("a b\nc\n")
    expected.write_text("c\nb a\n")
    names, species = ["a", "b", "c"], [0, 1, 2]
    if problem == "changed":
        expected.write_text("a\nb c\n")
    elif problem == "duplicate":
        partition.write_text("a b\nb c\n")
    elif problem == "missing":
        expected.write_text("a b\n")
    elif problem == "unknown":
        partition.write_text("a b\nc d\n")
    elif problem == "universe":
        names.append("a")
    elif problem == "overwrite":
        output.write_text("preserve")
    hits, src, dst, weights = ([0], [1], [4.]), [0], [1], [2.]
    def refine(clusters, *args, **kwargs):
        assert clusters == [[0, 1], [2]]
        assert args == ([], [], [], src, dst, species)
        assert kwargs == {"rbnh_scores": weights}
        return clusters
    def write(path, clusters, names):
        path.write_text("\n".join(" ".join(names[i] for i in group) for group in clusters) + "\n")
    replay = SimpleNamespace(read_index_clusters=lambda *args: [[0, 1], [2]],
        production_refinement_hits=lambda *args: ([], [], []), refine_cluster_indices=refine,
        write_clusters=write)
    args = (replay, read_partition, names, species, *hits, src, dst, weights, partition, expected, output)
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            module.reconstruct(*args)
        if problem == "overwrite":
            assert output.read_text() == "preserve"
    else:
        result = module.reconstruct(*args)
        assert result == dict(genes=3, groups=2, partition_equal=True, refinement_directed_hits=0)


def test_real_frozen_refinement_subprocess(tmp_path):
    root = Path(__file__).resolve().parents[2]
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if not launcher.is_dir():
        pytest.skip("Requires the retained publication launcher")
    (tmp_path / "input").write_text("a b\nc d\n")
    (tmp_path / "expected").write_text("d c\nb a\n")
    code = '''
import importlib.util, json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from benchmark_tools import replay_high_sensitivity as replay
from benchmark_tools.audit_historical_profile_ablation import read_partition
assert Path(replay.__file__).resolve().is_relative_to(Path(sys.argv[1]))
spec = importlib.util.spec_from_file_location('reconstruction', sys.argv[2])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
p = Path(sys.argv[3])
result = module.reconstruct(replay, read_partition, ['a','b','c','d'], [0,1,0,1],
    [], [], [], [0,2], [1,3], [2.,2.], p/'input', p/'expected', p/'output')
print(json.dumps(result))
'''
    done = subprocess.run([sys.executable, "-I", "-B", "-c", code, str(launcher),
        str(Path(module.__file__).resolve()), str(tmp_path)], cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert json.loads(done.stdout)["partition_equal"] is True


def test_cli_help(tmp_path):
    result = subprocess.run([sys.executable, "-I", str(Path(module.__file__).resolve()), "--help"],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
