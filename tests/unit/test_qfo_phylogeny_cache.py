import copy
import json

import pytest

from benchmark_tools.qfo_phylogeny_cache import cache_paths, copy_cache, record, select_cache


def fixture(tmp_path):
    source = tmp_path / "source"
    family = "Family0000001"
    for relative in cache_paths(family):
        path = source / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("content\n")
    raw = record(source / f"gene_trees/{family}.raw.nwk")
    checkpoint = dict(family_id=family, genes=["a", "b"], status="complete", schema_version=2,
                      tree_input_sha256="expected", raw_tree_sha256=raw["sha256"])
    (source / f"checkpoints/{family}.json").write_text(json.dumps(checkpoint))
    tool = tmp_path / "tool"
    tool.write_text("tool")
    records = [record(source / r) for r in cache_paths(family)]
    return source, {family: ("a", "b")}, records, lambda name, genes: "expected", [record(tool)]


def test_eligible_files_are_copied_without_trees_or_links(tmp_path):
    args = fixture(tmp_path)
    manifest = select_cache(*args)
    destination = tmp_path / "destination"
    copied = copy_cache(manifest, destination)
    assert len(copied) == 4
    assert not (destination / "species_tree_inference").exists()
    assert not list(destination.rglob("*.rooted.nwk"))
    raw = destination / "gene_trees/Family0000001.raw.nwk"
    original = args[0] / "gene_trees/Family0000001.raw.nwk"
    assert raw.stat().st_ino != original.stat().st_ino
    raw.write_text("changed")
    assert original.read_text() == "content\n"


@pytest.mark.parametrize("change", ["membership", "absent", "sequence_or_config"])
def test_changed_family_not_seeded(tmp_path, change):
    args = list(fixture(tmp_path))
    if change == "membership":
        args[1] = {"Family0000001": ("a", "c")}
    elif change == "absent":
        args[1] = {}
    else:
        args[3] = lambda name, genes: "different"
    result = select_cache(*args)
    assert result["included_families"] == []
    assert len(result["excluded_families"]) == 1
    assert result["files"] == []


def test_mutated_raw_or_tool_rejected(tmp_path):
    args = fixture(tmp_path)
    manifest = select_cache(*args)
    (args[0] / "gene_trees/Family0000001.raw.nwk").write_text("altered")
    with pytest.raises(ValueError):
        select_cache(*args)
    with pytest.raises(ValueError):
        copy_cache(manifest, tmp_path / "destination")


def test_changed_tool_rejected(tmp_path):
    args = fixture(tmp_path)
    (tmp_path / "tool").write_text("new tool")
    with pytest.raises(ValueError):
        select_cache(*args)


def test_unadmitted_alignment_rejected(tmp_path):
    args = list(fixture(tmp_path))
    args[2] = args[2][:-1]
    with pytest.raises(ValueError, match="Unadmitted"):
        select_cache(*args)


@pytest.mark.parametrize("name", ["../bad", "/absolute", "Family1", "Family0000001/other"])
def test_unsafe_family_id(name):
    with pytest.raises(ValueError):
        cache_paths(name)


def test_copy_refuses_extra_missing_duplicate_files_and_existing_destination(tmp_path):
    manifest = select_cache(*fixture(tmp_path))
    for files in (manifest["files"][:-1], manifest["files"] + [manifest["files"][0]],
                  [dict(relative="../bad", source=manifest["files"][0]["source"]) ]):
        changed = copy.deepcopy(manifest)
        changed["files"] = files
        with pytest.raises(ValueError):
            copy_cache(changed, tmp_path / "destination")
    destination = tmp_path / "destination"
    destination.mkdir()
    with pytest.raises(FileExistsError):
        copy_cache(manifest, destination)
