import os

import pytest

from benchmark_tools.isolated_numba_cache import fresh_cache


@pytest.mark.parametrize("prior", [None, "/existing/cache"])
@pytest.mark.parametrize("failure", [False, True])
def test_fresh_and_restore(tmp_path, monkeypatch, prior, failure):
    monkeypatch.setenv("NUMBA_CACHE_LOCATOR_CLASSES", "InTreeCacheLocator")
    if prior is None:
        monkeypatch.delenv("NUMBA_CACHE_DIR", raising=False)
    else:
        monkeypatch.setenv("NUMBA_CACHE_DIR", prior)
    cache = tmp_path / "cache"
    try:
        with fresh_cache(cache):
            assert cache.is_dir() and not list(cache.iterdir())
            assert os.environ["NUMBA_CACHE_DIR"] == str(cache)
            assert os.environ["NUMBA_CACHE_LOCATOR_CLASSES"] == "UserProvidedCacheLocator"
            (cache / "evidence").write_text("retained")
            if failure:
                raise RuntimeError("native failed")
    except RuntimeError:
        assert failure
    assert os.environ.get("NUMBA_CACHE_DIR") == prior
    assert os.environ["NUMBA_CACHE_LOCATOR_CLASSES"] == "InTreeCacheLocator"
    assert (cache / "evidence").read_text() == "retained"
    with pytest.raises(ValueError):
        with fresh_cache(cache):
            pass


def test_indirect_parent_rejected(tmp_path):
    (tmp_path / "real").mkdir()
    (tmp_path / "link").symlink_to(tmp_path / "real", target_is_directory=True)
    with pytest.raises(ValueError):
        with fresh_cache(tmp_path / "link/cache"):
            pass
