"""Opt-in pipeline metrics for reproducible performance benchmarks."""

from __future__ import annotations

import json
import os
import platform
import resource
import sys
import threading
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterator, Optional


def _process_rss_bytes(pid: int) -> int:
    """Return resident bytes for one Linux process, or zero if it exited."""
    try:
        fields = Path(f"/proc/{pid}/statm").read_text().split()
        return int(fields[1]) * os.sysconf("SC_PAGE_SIZE")
    except (FileNotFoundError, PermissionError, IndexError, ValueError, OSError):
        return 0


def _child_pids(pid: int) -> list[int]:
    try:
        text = Path(f"/proc/{pid}/task/{pid}/children").read_text().strip()
    except (FileNotFoundError, PermissionError, OSError):
        return []
    if not text:
        return []
    children = []
    for token in text.split():
        try:
            children.append(int(token))
        except ValueError:
            continue
    return children


def process_tree_rss_bytes(root_pid: int) -> int:
    """Return the current summed RSS of a process and all descendants."""
    total = 0
    pending = [root_pid]
    seen = set()
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += _process_rss_bytes(pid)
        pending.extend(_child_pids(pid))
    return total


class PipelineMetrics:
    """Collect stage timings and sampled process-tree memory when requested."""

    def __init__(self, output_path: Optional[str], sample_interval: float = 0.05):
        self.output_path = Path(output_path).resolve() if output_path else None
        self.sample_interval = sample_interval
        self.enabled = self.output_path is not None
        self.data: Dict[str, Any] = {}
        self._active_stage: Optional[str] = None
        self._stage_started = 0.0
        self._stage_cpu_started = (0.0, 0.0)
        self._stage_peaks: Dict[str, int] = {}
        self._overall_peak = 0
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None

    @staticmethod
    def _cpu_seconds() -> tuple[float, float]:
        own = resource.getrusage(resource.RUSAGE_SELF)
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        return own.ru_utime + children.ru_utime, own.ru_stime + children.ru_stime

    def __enter__(self) -> "PipelineMetrics":
        if not self.enabled:
            return self
        self.data = {
            "schema_version": 1,
            "status": "running",
            "command": [sys.executable, *sys.argv],
            "cwd": os.getcwd(),
            "pid": os.getpid(),
            "started_at_epoch_s": time.time(),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "stages": {},
            "counts": {},
        }
        self._sample()
        self._thread = threading.Thread(
            target=self._monitor,
            name="orthohmm-resource-monitor",
            daemon=True,
        )
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc_value, _traceback) -> bool:
        if not self.enabled:
            return False
        if self._active_stage is not None:
            self._finish_stage(self._active_stage)
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.sample_interval * 4))
        self._sample()
        user_s, system_s = self._cpu_seconds()
        completed = exc_type is None or (
            exc_type is SystemExit and getattr(exc_value, "code", 0) in (None, 0)
        )
        self.data.update({
            "status": "complete" if completed else "failed",
            "finished_at_epoch_s": time.time(),
            "wall_s": round(
                time.time() - float(self.data["started_at_epoch_s"]), 6
            ),
            "user_cpu_s": round(user_s, 6),
            "system_cpu_s": round(system_s, 6),
            "peak_process_tree_rss_bytes": self._overall_peak,
            "peak_process_tree_rss_gib": round(
                self._overall_peak / (1024 ** 3), 6
            ),
            "rss_measurement": (
                "sampled_sum_of_linux_proc_tree_rss"
                if Path("/proc/self/statm").exists()
                else "unavailable"
            ),
        })
        if exc_value is not None and not completed:
            self.data["error"] = f"{exc_type.__name__}: {exc_value}"
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.output_path.with_suffix(self.output_path.suffix + ".tmp")
        temporary.write_text(json.dumps(self.data, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, self.output_path)
        return False

    def _monitor(self) -> None:
        while not self._stop_event.wait(self.sample_interval):
            self._sample()

    def _sample(self) -> None:
        rss = process_tree_rss_bytes(os.getpid())
        self._overall_peak = max(self._overall_peak, rss)
        if self._active_stage is not None:
            self._stage_peaks[self._active_stage] = max(
                self._stage_peaks.get(self._active_stage, 0), rss
            )

    @contextmanager
    def stage(self, name: str) -> Iterator[None]:
        if not self.enabled:
            yield
            return
        if self._active_stage is not None:
            raise RuntimeError(
                f"metrics stage {self._active_stage!r} is already active"
            )
        self._active_stage = name
        self._stage_started = time.perf_counter()
        self._stage_cpu_started = self._cpu_seconds()
        self._sample()
        try:
            yield
        finally:
            self._finish_stage(name)

    def _finish_stage(self, name: str) -> None:
        self._sample()
        user_s, system_s = self._cpu_seconds()
        start_user_s, start_system_s = self._stage_cpu_started
        self.data["stages"][name] = {
            "wall_s": round(time.perf_counter() - self._stage_started, 6),
            "user_cpu_s": round(user_s - start_user_s, 6),
            "system_cpu_s": round(system_s - start_system_s, 6),
            "peak_process_tree_rss_bytes": self._stage_peaks.get(name, 0),
            "peak_process_tree_rss_gib": round(
                self._stage_peaks.get(name, 0) / (1024 ** 3), 6
            ),
        }
        self._active_stage = None

    def add_counts(self, **counts: Any) -> None:
        if self.enabled:
            self.data["counts"].update(counts)

    def add_metadata(self, **metadata: Any) -> None:
        if self.enabled:
            self.data.setdefault("metadata", {}).update(metadata)
