"""Serialize slow host observations off the native point-cadence thread."""

import math
import threading
import time


class PeriodicHostObserver:
    def __init__(self, monitor, period=30., join_timeout=15.):
        if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0
               for v in (period, join_timeout)):
            raise ValueError("Require positive finite observation and join bounds")
        self.monitor = monitor
        self.period = period
        self.join_timeout = join_timeout
        self.stop = threading.Event()
        self.thread = None
        self.error = None
        self.deadline = None

    def start(self, *, anchor=None):
        if self.thread is not None or self.stop.is_set():
            raise ValueError("Periodic observer cannot be restarted")
        now = time.monotonic()
        if anchor is not None and (type(anchor) not in (int, float)
                or not math.isfinite(anchor) or not 0 <= anchor <= now):
            raise ValueError("Require a finite observed monotonic anchor, not a future time")
        self.deadline = (now if anchor is None else anchor) + self.period
        self.thread = threading.Thread(target=self._observe, daemon=True)
        try:
            self.thread.start()
        except BaseException:
            self.stop.set()
            raise

    def _observe(self):
        try:
            deadline = self.deadline if self.deadline is not None else time.monotonic() + self.period
            while not self.stop.wait(max(0., deadline - time.monotonic())):
                self.monitor.observe()
                # A slow scan remains slow evidence, not a burst of catch-up scans.
                deadline += self.period
                now = time.monotonic()
                if deadline < now:
                    deadline = now + self.period
        except BaseException as error:
            self.error = error
            self.stop.set()

    def close(self):
        self.stop.set()
        if self.thread is not None and self.thread.ident is not None:
            self.thread.join(timeout=self.join_timeout)
            if self.thread.is_alive():
                raise TimeoutError("Owned host observation thread did not finish")
        if self.error is not None:
            raise RuntimeError("Owned host observation thread failed") from self.error

    def finish(self, launch, end):
        self.close()
        # All periodic scans finish before the post-command bracket and summary.
        self.monitor.observe()
        return self.monitor.summary(launch, end)
