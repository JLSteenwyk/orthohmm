import json

import pytest

from benchmark_tools import inspect_fastoma_runtime as runtime


def info():
    return {"ServerVersion": "28.2.2", "Driver": "overlay2", "CgroupDriver": "systemd",
            "CgroupVersion": "2", "DefaultRuntime": "runc", "KernelVersion": "kernel",
            "OperatingSystem": "Ubuntu", "OSType": "linux", "Architecture": "x86_64",
            "SecurityOptions": ["name=seccomp"]}


def install_probe(monkeypatch, data, context="default"):
    def run(argv):
        result = {"Client": {"Context": context}} if argv[1] == "version" else data
        return {"stdout": json.dumps(result)}
    monkeypatch.setattr(runtime, "run", run)


def test_snapshot_excludes_unrelated_docker_metadata(monkeypatch):
    data = {**info(), "ContainersRunning": 7, "HttpProxy": "not retained"}
    install_probe(monkeypatch, data)
    assert runtime.docker_state()["info"] == info()


@pytest.mark.parametrize("field,value", [("CgroupVersion", "1"), ("CgroupDriver", "cgroupfs"),
                                          ("DefaultRuntime", "nvidia"), ("OSType", "windows")])
def test_incompatible_configuration_rejected(monkeypatch, field, value):
    install_probe(monkeypatch, {**info(), field: value})
    with pytest.raises(ValueError, match="Unexpected Docker"):
        runtime.docker_state()


def test_nondefault_context_rejected(monkeypatch):
    install_probe(monkeypatch, info(), "remote")
    with pytest.raises(ValueError):
        runtime.docker_state()


@pytest.mark.parametrize("key", ["DOCKER_HOST", "DOCKER_CONTEXT"])
def test_remote_override_rejected_before_inspection(monkeypatch, key):
    monkeypatch.setenv(key, "remote")
    with pytest.raises(ValueError, match="local Docker"):
        runtime.inspect()


JAVA = ('openjdk version "17.0.18-internal" 2026-01-20\n'
        'OpenJDK Runtime Environment (build 17.0.18-internal+0-adhoc.conda.src)\n')


def test_exact_conda_java_build():
    runtime.validate_java(JAVA)


@pytest.mark.parametrize("text", ["", JAVA.replace("17.0.18", "17.0.19"),
                                  JAVA.replace("-internal", ""), JAVA.splitlines()[0],
                                  JAVA.replace("adhoc.conda.src", "different-build")])
def test_wrong_java_build_rejected(text):
    with pytest.raises(ValueError, match="Java"):
        runtime.validate_java(text)
