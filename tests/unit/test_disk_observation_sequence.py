import json
import tracemalloc

import pytest

from benchmark_tools.disk_observation_sequence import DiskObservations
from benchmark_tools.replay_threadripper_scaling import observations


def test_lazy_sequence_preserves_order_and_slices(tmp_path):
    points = DiskObservations(tmp_path)
    values = [dict(index=i) for i in range(5)]
    for point in values:
        points.append(point)
    assert list(points) == values
    assert points[-1] == values[-1]
    assert list(points[1:][::-1]) == values[1:][::-1]
    assert list(points[::2]) == values[::2]
    with pytest.raises(IndexError):
        points[5]
    with pytest.raises(ValueError):
        points[:].append({})
    with pytest.raises(FileExistsError):
        DiskObservations(tmp_path).append({})
    assert list(points) == values


def test_v3_records_bind_content_and_order(tmp_path):
    points = DiskObservations(tmp_path)
    points.append(dict(value=1))
    points.append(dict(value=2))
    paths = [points.path(i) for i in range(2)]
    report = dict(schema="threadripper_scaling_v3", point_records=points.records())
    assert list(observations(tmp_path, report, paths)) == list(points)
    with pytest.raises(ValueError):
        observations(tmp_path, report, paths[::-1])
    points.path(0).write_text(json.dumps(dict(value=99)))
    with pytest.raises(ValueError):
        observations(tmp_path, report, paths)


def test_legacy_embedded_points_still_verified(tmp_path):
    points = DiskObservations(tmp_path)
    points.append(dict(value=1))
    report = dict(schema="threadripper_scaling_v2", points=[dict(value=1)])
    assert observations(tmp_path, report, [points.path(0)]) == list(points)
    report["points"][0]["value"] = 2
    with pytest.raises(ValueError):
        observations(tmp_path, report, [points.path(0)])


def test_live_history_does_not_retain_payloads(tmp_path):
    points = DiskObservations(tmp_path)
    tracemalloc.start()
    try:
        for i in range(200):
            points.append(dict(index=i, payload="x" * 100_000))
        current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert len(points) == 200
    assert current < 1_000_000
    assert peak < 2_000_000
    assert points[-1]["index"] == 199


def test_reject_symlink(tmp_path):
    target = tmp_path / "target.json"
    target.write_text("{}")
    (tmp_path / "point_000000.json").symlink_to(target)
    with pytest.raises(ValueError):
        DiskObservations(tmp_path, 1)[0]
