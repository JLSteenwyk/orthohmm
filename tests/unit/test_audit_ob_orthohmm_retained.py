import pytest

from benchmark_tools.audit_ob_orthohmm_retained import resolve_records, timing


@pytest.mark.parametrize("paths", [["/absolute"], ["../escape"], ["a", "a"]])
def test_reject_manifest_paths(tmp_path, paths):
    with pytest.raises(ValueError):
        resolve_records(tmp_path, [dict(path=p, bytes=1, sha256="x") for p in paths])


def test_relative_paths_preserve_hashes(tmp_path):
    assert resolve_records(tmp_path, [dict(path="a", bytes=1, sha256="x")]) == [
        dict(path=str(tmp_path / "a"), bytes=1, sha256="x")]


def metrics():
    return dict(wall_s=10, user_cpu_s=20, system_cpu_s=1, peak_process_tree_rss_bytes=100,
                started_at_epoch_s=100, finished_at_epoch_s=110,
                rss_measurement="sampled_sum_of_linux_proc_tree_rss")


def test_time_semantics_preserved():
    assert timing(metrics())["wall_s"] == 10


@pytest.mark.parametrize("key,value", [("wall_s", 9), ("user_cpu_s", float("nan")),
                                     ("peak_process_tree_rss_bytes", -1), ("rss_measurement", "unique")])
def test_reject_time_mismatch(key, value):
    data = metrics()
    data[key] = value
    with pytest.raises(ValueError):
        timing(data)


def test_unavailable_zero_rss_is_not_admitted():
    data = dict(metrics(), rss_measurement="unavailable", peak_process_tree_rss_bytes=0)
    with pytest.raises(ValueError, match="memory semantics"):
        timing(data)
