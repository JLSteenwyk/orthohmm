"""Stop only a command launched with start_new_session=True on POSIX."""

import os
import signal
import time


def _signal_group(process, signum):
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        return False
    except PermissionError:
        # Darwin can report EPERM for a zombie-only group. Reap the leader,
        # but accept cleanup only if a fresh probe proves the group is gone.
        if process.poll() is not None:
            try:
                os.killpg(process.pid, 0)
            except ProcessLookupError:
                return False
        raise
    return True


def stop_owned_group(process, grace=5.):
    """Send TERM, then KILL after grace; leader exit alone is not group exit."""
    if not _signal_group(process, signal.SIGTERM):
        return process.wait()
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        if not _signal_group(process, 0):
            return process.wait()
        time.sleep(.02)
    _signal_group(process, signal.SIGKILL)
    return process.wait()
