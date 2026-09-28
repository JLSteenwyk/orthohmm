import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import run_cpm_readback_control as control


SOURCE = Path(__file__).resolve().parents[2] / "benchmark_tools/audit_historical_profile_ablation.py"


def test_reader_omits_module_imports(tmp_path):
    source = tmp_path / "source.py"
    source.write_text("raise RuntimeError('module must not execute')\n" +
                      "def read_partition(path, universe):\n    return universe\n")
    assert control.reader(source)(None, {"a"}) == {"a"}


@pytest.mark.parametrize("body", ["pass\n", "def read_partition(): pass\ndef read_partition(): pass\n",
                                  "@decorator\ndef read_partition(): pass\n"])
def test_reader_rejects_ambiguous_or_decorated_source(tmp_path, body):
    source = tmp_path / "source.py"
    source.write_text(body)
    with pytest.raises(ValueError):
        control.reader(source)


@pytest.mark.parametrize("partition", ["a a\nb\n", "a\na b\n", "a c\n", "a\n"])
def test_frozen_reader_rejects_invalid_membership(tmp_path, partition):
    path = tmp_path / "partition.txt"
    path.write_text(partition)
    with pytest.raises(ValueError):
        control.reader(SOURCE)(path, {"a", "b"})


def test_check_rejects_changed_bytes(tmp_path):
    path = tmp_path / "input"
    path.write_text("a")
    item = control.record(path)
    path.write_text("b")
    with pytest.raises(ValueError):
        control.check([item])


def test_real_no_site_worker(tmp_path):
    names = tmp_path / "names"
    initial = tmp_path / "initial"
    refined = tmp_path / "refined"
    names.write_text("a\nb\n")
    initial.write_text("a b\n")
    refined.write_text("a\nb\n")
    config = dict(checked_records=[control.record(p) for p in (SOURCE, names, initial, refined)],
                  source=str(SOURCE), names=str(names), initial=str(initial), refined=str(refined),
                  genes=2, groups=2)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    result = subprocess.run([sys.executable, "-B", "-S", control.__file__, "--worker", str(path)],
                            capture_output=True, text=True, timeout=30, check=True)
    child = json.loads(result.stdout)
    assert child["genes"] == child["groups"] == 2
    assert child["gc_enabled"] and not child["scientific_modules"]
    assert child["accuracy_admitted"] is False
