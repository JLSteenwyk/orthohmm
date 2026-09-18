import pytest

from benchmark_tools import probe_bracketed_cpu as module


def test_unscheduled_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(KeyError):
        module.run(tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("seconds", [0, -1, 4])
def test_unfrozen_cpu_duration_rejected(seconds):
    with pytest.raises(ValueError):
        module.burn(seconds)


def test_evaluation_checks_durations_before_screen():
    with pytest.raises(ValueError, match="native"):
        module.evaluate({}, {}, {"load": {"cpu_s": 1.9}}, None, 123, 100)
