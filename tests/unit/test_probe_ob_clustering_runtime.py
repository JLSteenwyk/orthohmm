import pytest

from benchmark_tools.probe_ob_clustering_runtime import verify_sources


def sources(prefix):
    return [{"path": f"/{prefix}/{name}", "sha256": name}
            for name in ("externals.py", "helpers.py", "leiden_worker.py")]


def test_equal_scientific_sources_allow_different_install_path():
    verify_sources(sources("old"), sources("new"))


def test_changed_scientific_source_refused():
    observed = sources("old")
    observed[0]["sha256"] = "changed"
    with pytest.raises(ValueError):
        verify_sources(observed, sources("new"))


def test_missing_source_refused():
    with pytest.raises(ValueError):
        verify_sources(sources("old")[:-1], sources("new"))
