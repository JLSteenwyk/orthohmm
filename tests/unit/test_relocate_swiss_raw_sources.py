import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.relocate_swiss_raw_sources import restore_inputs, stage, validate_records
from benchmark_tools import relocate_swiss_raw_sources as relocation


@pytest.fixture
def staged(tmp_path):
    original = tmp_path / "original"
    original.mkdir()
    a, b = original / "a", original / "b"
    a.write_bytes(b"first raw input")
    b.write_bytes(b"second raw input")
    items = [record(a), record(b), record(a)]
    source = original / "source.json"
    source.write_text(json.dumps({"records": items}))
    source_record = record(source)
    binding = stage(source, source_record["sha256"], "records", tmp_path / "staged")
    return items, source_record, Path(binding["path"]), binding["sha256"], original


def test_moved_tree_preserves_occurrences_without_originals(staged, tmp_path):
    items, source, binding, digest, original = staged
    moved = tmp_path / "moved"
    shutil.copytree(binding.parent, moved)
    shutil.rmtree(original)
    shutil.rmtree(binding.parent)
    restored, provenance = restore_inputs(items, source, moved / "bindings.json", digest)
    assert len(restored) == 3
    assert restored[0] == restored[2]
    assert len(list((moved / "inputs").iterdir())) == 2
    for old, new in zip(items, restored):
        assert {k: v for k, v in new.items() if k != "path"} == {
            k: v for k, v in old.items() if k != "path"}
        assert Path(new["path"]).parent == moved / "inputs"
    assert [e["original"] for e in provenance["entries"]] == items
    assert provenance["redistribution_authorized"] is False


def test_default_never_relocates_or_reads_unavailable_sources(staged):
    items, source, _, _, original = staged
    shutil.rmtree(original)
    restored, provenance = restore_inputs(items, source)
    assert restored is items
    assert provenance is None


@pytest.mark.parametrize("partial", ["path", "digest"])
def test_both_binding_arguments_required(staged, partial):
    items, source, binding, digest, _ = staged
    with pytest.raises(ValueError, match="both"):
        restore_inputs(items, source, binding if partial == "path" else None,
                       digest if partial == "digest" else None)


def test_independent_binding_digest_required(staged):
    items, source, binding, _, _ = staged
    with pytest.raises(ValueError, match="manifest changed"):
        restore_inputs(items, source, binding, "0" * 64)


@pytest.mark.parametrize("fault", ["omitted", "extra", "reordered", "sha", "size", "path",
    "artifact", "traversal", "absolute", "rights", "schema", "bool_schema", "source_sha",
    "source_size", "extra_field", "entry_field", "non_object"])
def test_rehashed_binding_still_cannot_change_frozen_inventory(staged, fault):
    items, source, binding, _, _ = staged
    data = json.loads(binding.read_text())
    if fault == "omitted":
        data["entries"].pop()
    elif fault == "extra":
        data["entries"].append(data["entries"][0])
    elif fault == "reordered":
        data["entries"].reverse()
        data["entries"][0], data["entries"][1] = data["entries"][1], data["entries"][0]
    elif fault in ("sha", "size", "path"):
        key, value = {"sha": ("sha256", "0" * 64), "size": ("bytes", 999),
                      "path": ("path", "/different/original")}[fault]
        data["entries"][0]["original"][key] = value
    elif fault in ("artifact", "traversal", "absolute"):
        data["entries"][0]["artifact"] = {"artifact": "inputs/wrong",
            "traversal": "../outside", "absolute": "/outside"}[fault]
    elif fault == "rights":
        data["redistribution_authorized"] = True
    elif fault in ("schema", "bool_schema"):
        data["schema_version"] = 2 if fault == "schema" else True
    elif fault in ("source_sha", "source_size"):
        data["source"]["sha256" if fault == "source_sha" else "bytes"] = "0" * 64 if fault == "source_sha" else 0
    elif fault == "extra_field":
        data["extra"] = 1
    elif fault == "entry_field":
        data["entries"][0]["extra"] = 1
    else:
        data = []
    binding.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        restore_inputs(items, source, binding, record(binding)["sha256"])


@pytest.mark.parametrize("fault", ["same_size", "truncated", "missing", "symlink", "parent_symlink"])
def test_actual_raw_bytes_and_paths_checked(staged, tmp_path, fault):
    items, source, binding, digest, _ = staged
    path = binding.parent / "inputs" / items[0]["sha256"]
    if fault == "same_size":
        path.write_bytes(b"X" + path.read_bytes()[1:])
    elif fault == "truncated":
        path.write_bytes(path.read_bytes()[:-1])
    elif fault == "missing":
        path.unlink()
    elif fault == "symlink":
        copy = tmp_path / "copy"
        shutil.copyfile(path, copy)
        path.unlink()
        path.symlink_to(copy)
    else:
        directory = binding.parent / "inputs"
        directory.rename(tmp_path / "linked_inputs")
        directory.symlink_to(tmp_path / "linked_inputs", target_is_directory=True)
    with pytest.raises(ValueError):
        restore_inputs(items, source, binding, digest)


def test_no_overwrite_even_if_destination_empty(staged, tmp_path):
    _, source, binding, _, _ = staged
    with pytest.raises(FileExistsError):
        stage(source["path"], source["sha256"], "records", binding.parent)
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileExistsError):
        stage(source["path"], source["sha256"], "records", empty)


def test_destination_symlink_parent_rejected(staged, tmp_path):
    _, source, _, _, _ = staged
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink parents"):
        stage(source["path"], source["sha256"], "records", link / "new")
    assert not (real / "new").exists()


def test_stage_retains_partial_attempt_when_original_changes_during_copy(staged, tmp_path, monkeypatch):
    items, source, _, _, _ = staged
    original_copy = shutil.copyfileobj
    def mutate_after_copy(src, dst, **kwargs):
        original_copy(src, dst, **kwargs)
        Path(items[0]["path"]).write_bytes(b"X" + Path(items[0]["path"]).read_bytes()[1:])
    monkeypatch.setattr(relocation.shutil, "copyfileobj", mutate_after_copy)
    output = tmp_path / "partial"
    with pytest.raises(ValueError, match="identity changed"):
        stage(source["path"], source["sha256"], "records", output)
    assert (output / "inputs").is_dir()
    assert not (output / "bindings.json").exists()


def test_binding_mutation_during_validation_rejected(staged, monkeypatch):
    items, source, binding, digest, _ = staged
    original_check = relocation.check
    def mutate_binding(item):
        original_check(item)
        binding.write_bytes(binding.read_bytes() + b" ")
    monkeypatch.setattr(relocation, "check", mutate_binding)
    with pytest.raises(ValueError, match="during validation"):
        restore_inputs(items, source, binding, digest)


@pytest.mark.parametrize("fault", ["wrong_source_hash", "changed_raw", "missing_raw", "key"])
def test_stage_admits_sources_before_creating_destination(staged, tmp_path, fault):
    items, source, _, _, _ = staged
    if fault == "changed_raw":
        Path(items[0]["path"]).write_bytes(b"changed")
    elif fault == "missing_raw":
        Path(items[0]["path"]).unlink()
    output = tmp_path / "rejected"
    with pytest.raises((ValueError, FileNotFoundError)):
        stage(source["path"], "0" * 64 if fault == "wrong_source_hash" else source["sha256"],
              "unexpected" if fault == "key" else "records", output)
    assert not output.exists()


@pytest.mark.parametrize("fault", ["empty", "negative", "bool_size", "short_sha", "uppercase_sha",
                                  "relative", "conflict", "extra_field"])
def test_malformed_source_records_rejected(staged, fault):
    items = json.loads(json.dumps(staged[0]))
    if fault == "empty":
        items = []
    elif fault in ("negative", "bool_size"):
        items[0]["bytes"] = -1 if fault == "negative" else True
    elif fault in ("short_sha", "uppercase_sha"):
        items[0]["sha256"] = "abc" if fault == "short_sha" else "A" * 64
    elif fault == "relative":
        items[0]["path"] = "relative"
    elif fault == "conflict":
        items[-1]["sha256"] = "0" * 64
    else:
        items[0]["extra"] = 1
    with pytest.raises(ValueError):
        validate_records(items)
