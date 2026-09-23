"""Exercise spawn's replay of the actual CPM entry point with competing trees."""

from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize("entrypoint,legacy", [("run_qfo_cpm_variant.py", False),
    ("run_qfo_cpm_control.py", False), ("run_qfo_cpm_variant.py", True)])
def test_spawn_preserves_native_package_over_executor(tmp_path, entrypoint, legacy):
    launcher = tmp_path / "native"
    executor = tmp_path / "executor"
    for tree in (launcher, executor):
        package = tree / "orthohmm"
        package.mkdir(parents=True)
        (package / "__init__.py").write_text("")
    scripts = executor / "benchmark_tools"
    scripts.mkdir()
    source = Path(__file__).resolve().parents[2] / "benchmark_tools" / entrypoint
    script = scripts / entrypoint
    script.write_bytes(source.read_bytes())
    if legacy:
        script.write_text(script.read_text().replace(
            'if __name__ != "__mp_main__":\n    sys.path.insert', 'if True:\n    sys.path.insert'))
    code = """
import __main__
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import sys
__main__.__file__ = sys.argv[1]
sys.path[:0] = [sys.argv[2], sys.argv[3]]
with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context('spawn')) as pool:
    actual = pool.submit(eval, "__import__('orthohmm').__file__").result(timeout=20)
assert actual == sys.argv[4], actual
print(actual)
"""
    expected = launcher / "orthohmm/__init__.py"
    result = subprocess.run([sys.executable, "-c", code, str(script), str(launcher), str(executor), str(expected)],
                            text=True, capture_output=True, timeout=30)
    if legacy:
        assert result.returncode != 0
        assert str(executor / "orthohmm/__init__.py") in result.stderr
    else:
        assert result.returncode == 0, result.stdout + result.stderr
        assert result.stdout.strip() == str(expected)
