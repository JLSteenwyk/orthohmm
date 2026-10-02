import os
from pathlib import Path
import shutil
import subprocess

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "benchmark_tools/install_ci_test_compiler.sh"


@pytest.fixture
def environment(tmp_path):
    prefix = tmp_path / "Homebrew $prefix's files"
    (prefix / "bin").mkdir(parents=True)
    commands = tmp_path / "commands"
    commands.mkdir()
    brew = commands / "brew"
    brew.write_text("#!/bin/sh\ncase \"$1 $2\" in\n"
                    "  'install gcc') exit 0 ;;\n"
                    "  '--prefix gcc') printf '%s\\n' \"$FAKE_BREW_PREFIX\" ;;\n"
                    "  *) exit 9 ;;\nesac\n")
    brew.chmod(0o700)
    temporary = tmp_path / "CI $temporary's files"
    temporary.mkdir()
    path = tmp_path / "github-path"
    env = dict(os.environ, PATH=str(commands) + os.pathsep + os.environ["PATH"],
               FAKE_BREW_PREFIX=str(prefix), RUNNER_TEMP=str(temporary), GITHUB_PATH=str(path))
    return prefix, temporary, path, env


def execute(env):
    return subprocess.run([shutil.which("bash"), str(SCRIPT)], env=env,
                          text=True, capture_output=True, timeout=10)


def compiler(prefix, version, *, executable=True):
    path = prefix / "bin" / ("gcc-" + version)
    path.write_text("#!/bin/sh\nprintf '%s\\n' 'fixture compiler version'\n")
    if executable:
        path.chmod(0o700)
    return path


def test_selects_numbered_driver_and_quotes_special_paths(environment):
    prefix, temporary, path, env = environment
    selected = compiler(prefix, "16")
    run = execute(env)
    assert run.returncode == 0, run.stderr
    tools = temporary / "orthohmm-ci-tools"
    assert (tools / "gcc").is_symlink() and (tools / "gcc").resolve() == selected
    assert path.read_text() == str(tools) + "\n"
    assert "fixture compiler version" in run.stdout


@pytest.mark.parametrize("case", ["missing", "ambiguous", "nonexecutable", "reused_tools", "version_failure"])
def test_invalid_driver_or_reused_tools_never_amends_path(environment, case):
    prefix, temporary, path, env = environment
    if case != "missing":
        compiler(prefix, "16", executable=case != "nonexecutable")
    if case == "ambiguous":
        compiler(prefix, "15")
    if case == "reused_tools":
        (temporary / "orthohmm-ci-tools").mkdir()
    if case == "version_failure":
        (prefix / "bin/gcc-16").write_text("#!/bin/sh\nexit 7\n")
    assert execute(env).returncode != 0
    assert not path.exists()


@pytest.mark.parametrize("name", ["RUNNER_TEMP", "GITHUB_PATH"])
def test_missing_ci_binding_fails_before_path_change(environment, name):
    _, _, path, env = environment
    env.pop(name)
    assert execute(env).returncode != 0
    assert not path.exists()
