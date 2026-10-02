import signal

import pytest

from benchmark_tools import owned_process_group as module


@pytest.mark.parametrize("stage", [signal.SIGTERM, 0, signal.SIGKILL])
@pytest.mark.parametrize("outcome", ["gone", "denied", "group_alive", "leader_alive"])
def test_permission_error_can_precede_waitable_exit(monkeypatch, stage, outcome):
    clock = [0.]
    signals, sleeps, injected = [], [], []
    failure = PermissionError("owned group permission denied")

    class Command:
        pid = 999999
        polls = 0
        waits = 0

        def poll(self):
            self.polls += 1
            return None if self.polls == 1 or outcome == "leader_alive" else -signal.SIGTERM

        def wait(self):
            self.waits += 1
            return -signal.SIGTERM

    def killpg(pid, signum):
        assert pid == command.pid
        signals.append(signum)
        if not injected and signum == stage:
            injected.append(signum)
            raise failure
        if injected:
            assert command.polls >= 2
            if outcome == "gone":
                raise ProcessLookupError("reaped group disappeared")
            if outcome == "group_alive":
                return
            raise failure

    def monotonic():
        clock[0] += .01
        return clock[0]

    command = Command()
    monkeypatch.setattr(module.os, "killpg", killpg)
    monkeypatch.setattr(module.time, "monotonic", monotonic)
    monkeypatch.setattr(module.time, "sleep", sleeps.append)
    grace = .05 if stage == 0 else 0.
    if outcome == "gone":
        assert module.stop_owned_group(command, grace) == -signal.SIGTERM
        assert signals[-2:] == [stage, 0]
        assert command.waits == 1
    else:
        with pytest.raises(PermissionError) as caught:
            module.stop_owned_group(command, grace)
        assert caught.value is failure
        assert command.waits == 0
        assert signals[-1] == (stage if outcome == "leader_alive" else 0)
    assert command.polls >= 2
    assert 0 < sum(sleeps) <= .1


@pytest.mark.parametrize("stage", [signal.SIGTERM, 0, signal.SIGKILL])
def test_waitability_near_bound_is_not_an_unbounded_wait(monkeypatch, stage):
    clock = [0.]
    injected = []
    polls = []

    class Command:
        pid = 999999

        def poll(self):
            polls.append(clock[0])
            return -signal.SIGTERM if clock[0] >= .07 else None

        def wait(self):
            assert len(polls) > 2
            return -signal.SIGTERM

    def killpg(pid, signum):
        assert pid == Command.pid
        if not injected and signum == stage:
            injected.append(signum)
            raise PermissionError("exit is not waitable yet")
        if injected:
            assert polls[-1] >= .07
            raise ProcessLookupError("gone")

    def sleep(seconds):
        clock[0] += seconds

    monkeypatch.setattr(module.os, "killpg", killpg)
    monkeypatch.setattr(module.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(module.time, "sleep", sleep)
    assert module.stop_owned_group(Command(), .01 if stage == 0 else 0.) == -signal.SIGTERM
    assert .07 <= clock[0] <= .11
