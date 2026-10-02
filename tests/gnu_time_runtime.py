"""Explicit GNU Time binding for smoke tests, never frozen benchmark execution."""

import os
from pathlib import Path
import shutil
import subprocess


def select_gnu_time():
    configured = os.environ.get("ORTHOHMM_TEST_GNU_TIME")
    candidates = ([configured] if configured is not None else
                  [shutil.which("gtime"), shutil.which("time"), "/usr/bin/time"])
    for candidate in dict.fromkeys(candidates):
        if not candidate:
            continue
        path = Path(candidate)
        if not path.is_absolute() or not path.is_file() or not os.access(path, os.X_OK):
            continue
        try:
            version = subprocess.run([str(path), "--version"], capture_output=True,
                                     text=True, check=False, timeout=5)
        except (OSError, subprocess.TimeoutExpired):
            continue
        if version.returncode == 0 and "GNU Time" in version.stdout:
            return str(path)
    if configured is not None:
        raise ValueError("ORTHOHMM_TEST_GNU_TIME must name an executable absolute GNU Time path")
    raise FileNotFoundError("GNU Time is required for native timing-wrapper smoke tests")


class GnuTimeSubprocess:
    """Replace only a module's literal timing executor, not global subprocess."""

    def __init__(self, executable, subprocess_module):
        self.executable = executable
        self.subprocess_module = subprocess_module

    def __getattr__(self, name):
        return getattr(self.subprocess_module, name)

    def _argv(self, args):
        if isinstance(args, (list, tuple)) and args and args[0] == "/usr/bin/time":
            return [self.executable, *args[1:]]
        return args

    def run(self, args, *positional, **keywords):
        return self.subprocess_module.run(self._argv(args), *positional, **keywords)

    def Popen(self, args, *positional, **keywords):
        return self.subprocess_module.Popen(self._argv(args), *positional, **keywords)
