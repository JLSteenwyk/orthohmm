import json

import pytest

from benchmark_tools.materialize_scaling_transfer import materialize, record


def fixture(tmp_path):
    bundle = tmp_path / "bundle"
    (bundle / "inputs").mkdir(parents=True)
    files = []
    for i in range(12):
        path = bundle / "inputs" / f"p{i}.fa"
        path.write_text(f">s{i}\nACDE\n")
        files.append({**record(path), "path": f"inputs/p{i}.fa"})
    manifest = {"inputs": files, "planned_runs": [], "datasets": [
        {"proteomes": n, "inputs": [f["path"] for f in files[:n]]} for n in (4, 8, 12)]}
    path = bundle / "manifest.json"
    path.write_text(json.dumps(manifest))
    return bundle, record(path)["sha256"]


def test_materialized_inputs_match(tmp_path):
    bundle, digest = fixture(tmp_path)
    result = materialize(bundle, tmp_path / "output", digest)
    assert [len(d["inputs"]) for d in result["datasets"]] == [4, 8, 12]
    assert result["inference_started"] is False
    assert (tmp_path / "output/12/p0.fa").read_bytes() == (bundle / "inputs/p0.fa").read_bytes()
    with pytest.raises(FileExistsError):
        materialize(bundle, tmp_path / "output", digest)


@pytest.mark.parametrize("change", ["hash", "extra", "symlink", "manifest"])
def test_reject_changed_transfer_before_copy(tmp_path, change):
    bundle, digest = fixture(tmp_path)
    path = bundle / "inputs/p0.fa"
    if change == "hash":
        path.write_text("changed")
    elif change == "extra":
        (bundle / "inputs/extra.fa").write_text("extra")
    elif change == "symlink":
        path.unlink()
        path.symlink_to(bundle / "inputs/p1.fa")
    else:
        digest = "wrong"
    with pytest.raises(ValueError):
        materialize(bundle, tmp_path / "output", digest)
    assert not (tmp_path / "output").exists()
