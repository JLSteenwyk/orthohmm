import pytest

from benchmark_tools import probe_scaling_collector as module


@pytest.mark.parametrize("fault", [None, "exit", "timeout", "missing", "periodic", "bracket"])
def test_probe_requires_native_success_and_complete_observation(tmp_path, monkeypatch, fault):
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    host = dict(successful_snapshots=3, command_bracketed_by_samples=True,
                controlled_workload_verified=False)
    native = dict(exit_code=0, timed_out=False)
    if fault == "exit": native["exit_code"] = 7
    if fault == "timeout": native["timed_out"] = True
    if fault == "missing": host = None
    if fault == "periodic": host["successful_snapshots"] = 2
    if fault == "bracket": host["command_bracketed_by_samples"] = False
    def measure(command, path, job, cpus, memory, timeout, interval):
        assert command == ["/bin/sleep", "35"] and job == 123
        assert (cpus, memory, timeout, interval) == (20, 96 * 1024**3, 85800, 1.)
        return dict(native=native, native_wall_s=35.)
    monkeypatch.setattr(module, "measure", measure)
    monkeypatch.setattr(module, "replay", lambda *a: dict(host_process_replay=host))
    output = tmp_path / "probe"
    if fault:
        with pytest.raises(ValueError): module.run(output)
        assert not (output / "probe.json").exists()
    else:
        result = module.run(output)
        assert result["status"] == "collector_probe_completed"
        assert result["scientific_timings_admitted"] is False
    assert (output / "replay.json").exists()
