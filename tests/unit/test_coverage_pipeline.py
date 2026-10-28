"""Exercise coverage collection and combination, not just config spelling."""

import configparser
from importlib import metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from coverage import CoverageData
from packaging.markers import default_environment
from packaging.requirements import Requirement
import pytest


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("version,expected", [("3.10", "7.16.2"), ("3.11", "7.17.0"),
                                            ("3.12", "7.17.0"), ("3.13", "7.17.0")])
def test_coverage_pins_select_one_supported_version(version, expected):
    environment = dict(default_environment(), python_version=version, python_full_version=version + ".0")
    requirements = [Requirement(line) for line in (ROOT / "tests/requirements.txt").read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith("#")]
    selected = [r for r in requirements if r.name == "coverage"
                and (r.marker is None or r.marker.evaluate(environment))]
    assert len(selected) == 1 and str(selected[0].specifier) == "==" + expected


def test_current_interpreter_uses_its_declared_coverage_pin():
    requirements = [Requirement(line) for line in (ROOT / "tests/requirements.txt").read_text().splitlines()
                    if line.strip() and not line.lstrip().startswith("#")]
    selected = [r for r in requirements if r.name == "coverage"
                and (r.marker is None or r.marker.evaluate())]
    assert len(selected) == 1 and selected[0].specifier.contains(metadata.version("coverage"))


def clean_environment():
    return {key: value for key, value in os.environ.items()
            if not key.startswith(("COVERAGE_", "COV_CORE_"))}


def checked(command, directory, environment):
    result = subprocess.run(command, cwd=directory, env=environment,
                            capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def test_current_config_keeps_product_and_benchmark_sources_visible():
    config = configparser.ConfigParser()
    config.read(ROOT / ".coveragerc")
    assert config.get("run", "source").strip() == "."
    assert set(config.get("run", "patch").split()) == {"subprocess", "_exit"}
    assert config.getboolean("run", "relative_files")
    omissions = config.get("run", "omit").split()
    assert omissions == ["tests/*", "setup.py", "orthohmm-runner.py"]
    assert "orthohmm/__main__.py" not in (ROOT / "codecov.yml").read_text()
    workflow = (ROOT / ".github/workflows/ci.yml").read_text()
    assert "coverage.json" in workflow and "coverage.xml" in workflow
    assert "NUMBA_JIT_COVERAGE" not in workflow


@pytest.mark.parametrize("flush_abrupt_exit", [False, True])
def test_subprocess_coverage_survives_changed_cwd_and_abrupt_exit(tmp_path, flush_abrupt_exit):
    app = tmp_path / "app"
    elsewhere = tmp_path / "different working directory"
    app.mkdir()
    elsewhere.mkdir()
    worker = app / "worker.py"
    worker.write_text("import os, sys\n"
                      "def run():\n"
                      "    marker = 7\n"
                      "    if sys.argv[1] == 'abrupt':\n"
                      "        os._exit(0)\n"
                      "    return marker\n"
                      "run()\n")
    driver = app / "driver.py"
    driver.write_text("import subprocess, sys\n"
                      "for mode in ['normal', 'abrupt']:\n"
                      "    subprocess.run([sys.executable, sys.argv[1], mode], cwd=sys.argv[2], check=True)\n")
    config = configparser.ConfigParser()
    config.read(ROOT / ".coveragerc")
    config.set("run", "source", str(app))
    config.set("run", "patch", "subprocess\n_exit" if flush_abrupt_exit else "subprocess")
    config_path = tmp_path / "coverage.ini"
    with config_path.open("w") as handle:
        config.write(handle)
    data_path = tmp_path / ".coverage"
    environment = clean_environment()
    environment["COVERAGE_FILE"] = str(data_path)
    command = [sys.executable, "-m", "coverage"]
    checked([*command, "run", "--rcfile", str(config_path), str(driver), str(worker), str(elsewhere)],
            tmp_path, environment)
    checked([*command, "combine", "--rcfile", str(config_path)], tmp_path, environment)
    data = CoverageData(basename=str(data_path))
    data.read()
    filename = next(name for name in data.measured_files() if Path(name).name == "worker.py")
    lines = data.lines(filename)
    assert 3 in lines and 6 in lines
    assert (5 in lines) is flush_abrupt_exit


@pytest.mark.parametrize("missing_application_line", [False, True])
def test_make_combines_unit_and_integration_without_discarding_either(tmp_path, missing_application_line):
    if shutil.which("make") is None:
        pytest.skip("GNU/BSD make required")
    app = tmp_path / "orthohmm"
    app.mkdir()
    source = app / "tiny.py"
    source.write_text("first = 1\nsecond = 2\n" + ("untested = 3\n" if missing_application_line else ""))
    shutil.copyfile(ROOT / "Makefile", tmp_path / "Makefile")
    shutil.copyfile(ROOT / ".coveragerc", tmp_path / ".coveragerc")
    for label, line in (("unit", 1), ("integration", 2)):
        data = CoverageData(basename=str(tmp_path / (".coverage-" + label)))
        data.add_lines({"orthohmm/tiny.py": {line}})
        data.write()
    environment = clean_environment()
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", "")
    result = subprocess.run(["make", "coverage.combine"], cwd=tmp_path, env=environment,
                            capture_output=True, text=True, timeout=90)
    assert (result.returncode != 0) is missing_application_line, result.stdout + result.stderr
    if missing_application_line:
        assert "Coverage failure" in result.stdout
    combined = CoverageData(basename=str(tmp_path / ".coverage"))
    combined.read()
    assert set(combined.lines("orthohmm/tiny.py")) == {1, 2}
    for label, line in (("unit", 1), ("integration", 2)):
        original = CoverageData(basename=str(tmp_path / (".coverage-" + label)))
        original.read()
        assert original.lines("orthohmm/tiny.py") == [line]
    report = json.loads((tmp_path / "coverage.json").read_text())
    summary = report["files"]["orthohmm/tiny.py"]["summary"]
    assert summary["covered_lines"] == 2
    assert summary["num_statements"] == (3 if missing_application_line else 2)
    assert (tmp_path / "coverage.xml").is_file()
