"""Keep native temporary files under a fresh per-run persistent directory."""

from contextlib import contextmanager
import os
from pathlib import Path


@contextmanager
def fresh_tmp(directory):
    directory = Path(directory)
    if (not directory.is_absolute() or directory.resolve() != directory
            or directory.exists() or directory.is_symlink()
            or directory.is_relative_to("/dev/shm")):
        raise ValueError("Require a fresh direct persistent scratch path")
    directory.mkdir(parents=True, exist_ok=False)
    keys = ("TMPDIR", "TMP", "TEMP")
    previous = {key: os.environ.get(key) for key in keys}
    os.environ.update({key: str(directory) for key in keys})
    try:
        yield directory
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
