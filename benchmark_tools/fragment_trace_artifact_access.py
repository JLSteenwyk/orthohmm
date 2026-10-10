"""Read-only access to explicitly indexed copies of frozen logical artifacts."""

import hashlib
import os
from pathlib import Path, PurePosixPath


def identity(path):
    data = Path(path).read_bytes()
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


class ArtifactAccess:
    def __init__(self, root, rows):
        self.root = Path(root).resolve(strict=True)
        self.rows, self.reverse = {}, {}
        for row in rows:
            logical = row["logical_path"]
            target = PurePosixPath(row["target"])
            if (not PurePosixPath(logical).is_absolute() or str(PurePosixPath(logical)) != logical
                    or ".." in PurePosixPath(logical).parts or logical in self.rows
                    or target.is_absolute() or ".." in target.parts or str(target) != row["target"]):
                raise ValueError("Invalid or duplicate artifact map")
            physical = self.root / target
            if str(physical) in self.reverse:
                raise ValueError("Multiple logical paths share one target")
            self.rows[logical] = row
            self.reverse[str(physical)] = logical

    def path(self, value):
        return ArtifactPath(self, value)

    def physical(self, logical):
        if logical not in self.rows:
            raise ValueError("Unindexed artifact: " + logical)
        row = self.rows[logical]
        target = self.root / row["target"]
        if not target.is_file() or target.is_symlink() or not target.resolve().is_relative_to(self.root):
            raise ValueError("Missing or escaped copied artifact")
        if identity(target) != {k:row[k] for k in ("bytes", "sha256")}:
            raise ValueError("Copied artifact differs: " + logical)
        return target

    def verify(self):
        for logical in self.rows:
            self.physical(logical)


class ArtifactPath(os.PathLike):
    """Keep historical labels visible while reads use verified copied files."""

    def __init__(self, access, value):
        self.access = access
        value = str(value)
        self.logical = PurePosixPath(access.reverse.get(value, value))
        if not self.logical.is_absolute() or ".." in self.logical.parts:
            raise ValueError("Require absolute logical artifact path")

    def __str__(self):
        return str(self.logical)

    def __fspath__(self):
        # NumPy accepts path-like values and must never open the historical path.
        return str(self.access.physical(str(self)))

    def __truediv__(self, name):
        child = PurePosixPath(name)
        if child.is_absolute() or ".." in child.parts:
            raise ValueError("Unsafe logical child")
        return ArtifactPath(self.access, self.logical / child)

    @property
    def name(self):
        return self.logical.name

    @property
    def stem(self):
        return self.logical.stem

    def resolve(self):
        return self

    def with_name(self, name):
        if PurePosixPath(name).name != name or name in ("", ".", ".."):
            raise ValueError("Unsafe replacement name")
        return ArtifactPath(self.access, self.logical.with_name(name))

    def glob(self, pattern):
        if PurePosixPath(pattern).is_absolute() or ".." in PurePosixPath(pattern).parts:
            raise ValueError("Unsafe logical glob")
        prefix = str(self.logical) + "/"
        for name in sorted(self.access.rows):
            if not name.startswith(prefix):
                continue
            tail = PurePosixPath(name[len(prefix):])
            if "/" not in pattern and len(tail.parts) != 1:
                continue
            if tail.match(pattern):
                yield ArtifactPath(self.access, name)

    def open(self, mode="r", *args, **kwargs):
        if mode not in ("r", "rt", "rb"):
            raise ValueError("Artifact access is read-only")
        return self.access.physical(str(self)).open(mode, *args, **kwargs)

    def read_bytes(self):
        return self.access.physical(str(self)).read_bytes()

    def read_text(self, *args, **kwargs):
        return self.access.physical(str(self)).read_text(*args, **kwargs)
