import threading
from types import SimpleNamespace

import pytest

from benchmark_tools.periodic_host_observer import PeriodicHostObserver


@pytest.mark.parametrize("period,timeout", [(0, 1), (1, 0), (-1, 1), (float('inf'), 1), (1, float('nan')), (True, 1)])
def test_invalid_bounds(period, timeout):
    with pytest.raises(ValueError):
        PeriodicHostObserver(None, period, timeout)


def test_close_before_start_and_restart_refused():
    observer = PeriodicHostObserver(None)
    observer.close()
    with pytest.raises(ValueError, match="restarted"):
        observer.start()


def test_slow_scan_does_not_block_calling_thread_and_final_bracket_is_serial():
    entered, release = threading.Event(), threading.Event()
    calls, active = [], []
    def observe():
        assert not active
        active.append(True)
        calls.append(threading.get_ident())
        try:
            if len(calls) == 1:
                entered.set()
                assert release.wait(2)
        finally:
            active.pop()
    monitor = SimpleNamespace(observe=observe, summary=lambda a, b: dict(start=a, end=b, calls=len(calls)))
    observer = PeriodicHostObserver(monitor, .01, 2)
    observer.start()
    try:
        assert entered.wait(2)
        assert observer.thread.is_alive()
        assert calls == [observer.thread.ident]
        with pytest.raises(ValueError, match="restarted"):
            observer.start()
        # Stop before releasing the scan so there is no timing-based race for another one.
        observer.stop.set()
        release.set()
        result = observer.finish(1., 2.)
        assert result == dict(start=1., end=2., calls=2)
        assert calls[1] == threading.get_ident()
        assert not observer.thread.is_alive()
        observer.close()
    finally:
        release.set()
        observer.close()


def test_unexpected_worker_exception_is_not_quiet_evidence():
    entered = threading.Event()
    def fail():
        entered.set()
        raise ValueError("synthetic scan failure")
    observer = PeriodicHostObserver(SimpleNamespace(observe=fail), .01, 2)
    observer.start()
    assert entered.wait(2)
    with pytest.raises(RuntimeError, match="thread failed") as error:
        observer.close()
    assert isinstance(error.value.__cause__, ValueError)
    assert not observer.thread.is_alive()


def test_join_timeout_does_not_summarize_or_start_another_scan():
    entered, release = threading.Event(), threading.Event()
    calls = []
    def block():
        entered.set()
        assert release.wait(2)
    observer = PeriodicHostObserver(SimpleNamespace(observe=block, summary=lambda *a: calls.append('summary')),
                                    .01, .01)
    observer.start()
    try:
        assert entered.wait(2)
        with pytest.raises(TimeoutError, match="did not finish"):
            observer.finish(1., 2.)
        assert not calls
        assert observer.stop.is_set()
    finally:
        release.set()
        observer.thread.join(2)
        observer.close()


def test_final_scan_failure_propagates():
    observer = PeriodicHostObserver(SimpleNamespace(observe=lambda: (_ for _ in ()).throw(OSError('final'))))
    with pytest.raises(OSError, match="final"):
        observer.finish(1., 2.)


def test_failed_thread_start_can_be_closed(monkeypatch):
    observer = PeriodicHostObserver(None)
    monkeypatch.setattr(threading.Thread, "start", lambda self: (_ for _ in ()).throw(RuntimeError('start')))
    with pytest.raises(RuntimeError, match="start"):
        observer.start()
    observer.close()


def test_slow_scan_skips_missed_periods_without_catchup_burst(monkeypatch):
    from benchmark_tools import periodic_host_observer as module
    clock, waits, calls = [0.], [], []
    monkeypatch.setattr(module.time, "monotonic", lambda: clock[0])
    def observe():
        calls.append(True)
        clock[0] = 100.
    class Stop:
        def wait(self, seconds):
            waits.append(seconds)
            return len(waits) == 2
    observer = PeriodicHostObserver(SimpleNamespace(observe=observe), 30.)
    observer.stop = Stop()
    observer._observe()
    assert waits == [30., 30.]
    assert calls == [True]


@pytest.mark.parametrize('anchor', [True, -1., float('nan'), float('inf'), 101., '100'])
def test_invalid_or_future_anchor_refused_before_thread_start(monkeypatch, anchor):
    from benchmark_tools import periodic_host_observer as module
    monkeypatch.setattr(module.time, 'monotonic', lambda: 100.)
    observer = PeriodicHostObserver(None)
    with pytest.raises(ValueError, match='monotonic anchor'):
        observer.start(anchor=anchor)
    assert observer.thread is None and not observer.stop.is_set()


def test_release_delay_does_not_shift_first_snapshot_deadline(monkeypatch):
    from benchmark_tools import periodic_host_observer as module
    clock, waits, scans = [3.2], [], []
    monkeypatch.setattr(module.time, 'monotonic', lambda: clock[0])
    class Thread:
        ident = None
        def __init__(self, **kwargs):
            pass
        def start(self):
            pass
    monkeypatch.setattr(module.threading, 'Thread', Thread)
    observer = PeriodicHostObserver(SimpleNamespace(observe=lambda: scans.append(clock[0])), 30.)
    observer.start(anchor=0.)
    assert observer.deadline == 30.
    # Initial inventory plus environmental handoff delayed native release by 13.2s.
    clock[0] = 13.2
    class Stop:
        def wait(self, seconds):
            waits.append(seconds)
            clock[0] += seconds
            return len(waits) == 2
    observer.stop = Stop()
    observer._observe()
    assert waits == [16.8, 30.]
    assert scans == [30.]
