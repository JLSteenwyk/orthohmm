"""Fresh persistent JIT caches with restoration of the caller environment."""

from contextlib import contextmanager
import os
from pathlib import Path


@contextmanager
def fresh_cache(directory):
    directory = Path(directory)
    if (not directory.is_absolute() or directory.resolve() != directory
            or directory.exists() or directory.is_symlink()):
        raise ValueError("Require a fresh direct absolute Numba cache path")
    directory.mkdir(parents=True, exist_ok=False)
    overrides = {"NUMBA_CACHE_DIR": str(directory),
                 "NUMBA_CACHE_LOCATOR_CLASSES": "UserProvidedCacheLocator"}
    previous = {key: os.environ.get(key) for key in overrides}
    os.environ.update(overrides)
    try:
        yield directory
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
