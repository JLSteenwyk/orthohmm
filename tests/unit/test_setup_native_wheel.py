from pathlib import Path
import runpy

import setuptools


def test_ctypes_distribution_is_not_platform_independent(monkeypatch):
    captured = {}
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: captured.update(kwargs))
    runpy.run_path(str(Path(__file__).resolve().parents[2] / "setup.py"))
    distribution = captured["distclass"]()
    assert distribution.has_ext_modules()
    assert not distribution.is_pure()
    command = distribution.get_command_obj("bdist_wheel")
    command.ensure_finalized()
    assert command.root_is_pure is False
    assert command.get_tag()[2] != "any"
    assert captured["entry_points"]["console_scripts"] == [
        "orthohmm = orthohmm.orthohmm:main"
    ]
