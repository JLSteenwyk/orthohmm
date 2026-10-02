import os

import pytest

from tests.conftest import isolated_launcher_environment


@pytest.mark.parametrize("exceptional", [False, True])
@pytest.mark.parametrize("teardown_order", ["fixture_first", "monkeypatch_first"])
def test_restores_complete_environment(monkeypatch, exceptional, teardown_order):
    monkeypatch.setenv("PATH", "/synthetic/original/toolchain")
    monkeypatch.setenv("ORTHOHMM_ENV_KEEP", "original")
    before = dict(os.environ)
    patch = pytest.MonkeyPatch()
    fixture = isolated_launcher_environment.__wrapped__()
    next(fixture)
    try:
        patch.setenv("ORTHOHMM_ENV_MONKEYPATCH", "temporary")
        os.environ["PATH"] = "/synthetic/launcher/toolchain"
        os.environ.pop("ORTHOHMM_ENV_KEEP")
        os.environ["ORTHOHMM_ENV_ADDED"] = "temporary"
        if teardown_order == "monkeypatch_first":
            patch.undo()
        if exceptional:
            with pytest.raises(RuntimeError, match="synthetic failure"):
                fixture.throw(RuntimeError("synthetic failure"))
        else:
            fixture.close()
        assert dict(os.environ) == before
        patch.undo()
        assert dict(os.environ) == before
    finally:
        fixture.close()
        patch.undo()
