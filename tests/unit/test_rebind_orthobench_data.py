import json
from pathlib import Path

import pytest

from benchmark_tools import rebind_orthobench_data as module


def fixture(tmp_path):
    root = tmp_path / "acquired"
    data = dict(dataset="orthobench", genes=251378, upstream_revision="synthetic")
    for (role, directory), count in zip(module.ROLES.items(), (12, 70, 11)):
        folder = root / directory
        folder.mkdir(parents=True, exist_ok=True)
        data[role] = []
        for i in range(count):
            path = folder / f"{i}.txt"
            path.write_text(f"synthetic {role} {i}\n")
            row = module.record(path)
            row["path"] = str(Path("/unavailable-original") / directory / path.name)
            data[role].append(row)
    source = tmp_path / "manifest.json"
    source.write_text(json.dumps(data))
    return source, root, data


def test_rebinding_needs_no_original_paths_and_preserves_order(tmp_path):
    source, root, old = fixture(tmp_path)
    result = module.rebind(source, module.record(source)["sha256"], root, tmp_path / "out")
    new = json.loads(Path(result["manifest"]["path"]).read_text())
    assert result["files"] == 93
    assert result["inference_rerun"] is result["publication_ready"] is False
    for role in module.ROLES:
        assert [(r["bytes"], r["sha256"]) for r in new[role]] == [
            (r["bytes"], r["sha256"]) for r in old[role]]
        assert all(Path(r["path"]).is_relative_to(root) for r in new[role])
    assert new["upstream_revision"] == old["upstream_revision"]
    with pytest.raises(FileExistsError):
        module.rebind(source, result["original_manifest"]["sha256"], root, tmp_path / "out")


@pytest.mark.parametrize("defect", ["digest", "changed", "missing", "duplicate", "role", "escape", "root", "count", "symlink"])
def test_invalid_inputs_do_not_publish_manifest(tmp_path, defect):
    source, root, data = fixture(tmp_path)
    path = root / "BENCHMARKS/Input/0.txt"
    if defect == "changed":
        path.write_text("changed")
    elif defect == "missing":
        path.unlink()
    elif defect == "duplicate":
        data["fasta"][1] = data["fasta"][0]
    elif defect == "role":
        data["fasta"][0]["path"] = "/original/BENCHMARKS/RefOGs/0.txt"
    elif defect == "escape":
        data["fasta"][0]["path"] = "/original/../BENCHMARKS/Input/0.txt"
    elif defect == "root":
        data["fasta"][0]["path"] = "/different/BENCHMARKS/Input/0.txt"
    elif defect == "count":
        data["fasta"].pop()
    elif defect == "symlink":
        copy = tmp_path / "outside"
        path.rename(copy)
        path.symlink_to(copy)
    source.write_text(json.dumps(data))
    digest = "0" * 64 if defect == "digest" else module.record(source)["sha256"]
    with pytest.raises((ValueError, FileNotFoundError)):
        module.rebind(source, digest, root, tmp_path / "out")
    assert not (tmp_path / "out").exists()
