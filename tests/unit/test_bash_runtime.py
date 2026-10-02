import shutil
import subprocess

import pytest


def bash():
    shell = shutil.which("bash")
    assert shell is not None, "Batch handoff tests require Bash"
    version = subprocess.run([shell, "--version"], check=True, text=True, capture_output=True)
    return shell, version.stdout.splitlines()[0]


@pytest.mark.parametrize("guard", [
    "false",
    '[[ "wrong_commit" == "expected_commit" ]]',
    '[[ "HEAD" =~ ^[0-9a-f]{40}$ ]]',
    '[[ "pending" =~ ^[0-9]+$ ]]',
    'INDEX=1\n[[ "$INDEX" =~ ^[0246]$ ]]',
])
def test_batch_shell_stops_on_failed_guard(guard):
    shell, version = bash()
    result = subprocess.run([shell, "-c", "set -euo pipefail\n" + guard +
                             "\nprintf 'guard-bypassed\\n'\n"], text=True, capture_output=True)
    assert result.returncode == 1 and result.stdout == "", (
        f"Unsupported batch guard semantics: {shell}, {version}, "
        f"returncode={result.returncode}, stdout={result.stdout!r}, stderr={result.stderr!r}")


def test_batch_shell_accepts_valid_guard():
    shell, version = bash()
    result = subprocess.run([shell, "-c", 'set -euo pipefail\n[[ "123" =~ ^[0-9]+$ ]]\n'
                             "printf 'guard-accepted\\n'\n"], text=True, capture_output=True)
    assert result.returncode == 0 and result.stdout == "guard-accepted\n", (shell, version, result)
