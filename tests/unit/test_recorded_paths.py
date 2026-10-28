import copy
import hashlib
import json

from tests.recorded_paths import RECORDED_ROOT, load_json, relocate


def test_relocation_keeps_original_content_and_identity_fields(tmp_path):
    value = dict(path=RECORDED_ROOT + "/source.py", bytes=42, sha256="f" * 64,
                 nested=[RECORDED_ROOT, dict(other="/unrelated/input")], note="original prose")
    original = copy.deepcopy(value)
    result = relocate(value, tmp_path)
    assert result["path"] == str(tmp_path / "source.py")
    assert result["nested"] == [str(tmp_path), dict(other="/unrelated/input")]
    assert result["bytes"] == value["bytes"] and result["sha256"] == value["sha256"]
    assert value == original and result is not value


def test_only_exact_path_prefix_is_relocated(tmp_path):
    for text in (RECORDED_ROOT + "_other/file", "prefix " + RECORDED_ROOT, "relative/file"):
        assert relocate(text, tmp_path) == text


def test_loading_never_rewrites_retained_input(tmp_path):
    path = tmp_path / "input.json"
    path.write_text(json.dumps(dict(path=RECORDED_ROOT + "/input")))
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    assert load_json(path)["path"].endswith("/input")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
