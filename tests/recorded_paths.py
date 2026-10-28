"""Relocate paths in test copies without changing retained hashes or files."""

import json
import hashlib
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
RECORDED_ROOT = "/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm"


def relocate(value, root=ROOT):
    if isinstance(value, dict):
        return {key: relocate(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [relocate(item, root) for item in value]
    if isinstance(value, str):
        if value == RECORDED_ROOT:
            return str(root)
        if value.startswith(RECORDED_ROOT + "/"):
            return str(root) + value[len(RECORDED_ROOT):]
    return value


def load_json(path):
    return relocate(json.loads(Path(path).read_text()))


def loads(text):
    return relocate(json.loads(text))


def check_record(reference, *, revision=None):
    path = Path(relocate(reference["path"]))
    if not path.is_absolute():
        path = ROOT / path
    data = (path.read_bytes() if revision is None else subprocess.check_output(
        ["git", "show", revision + ":" + path.relative_to(ROOT).as_posix()], cwd=ROOT))
    assert (len(data), hashlib.sha256(data).hexdigest()) == (reference["bytes"], reference["sha256"])
    return path
