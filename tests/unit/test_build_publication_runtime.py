import json

import pytest

from benchmark_tools.build_publication_runtime import COMMIT, KERNELS, record, verify_runtime


def fixture_manifest(tmp_path):
    root = tmp_path / "root"
    csrc = root / "orthohmm/search/csrc"
    csrc.mkdir(parents=True)
    (root / "setup.py").write_text("frozen setup")
    for name in KERNELS:
        (csrc / (name + ".c")).write_text("source")
        (csrc / (name + ".so")).write_bytes(b"binary")
    data = {"status": "complete", "commit": COMMIT, "root": str(root),
            "profile_probe": {"status": "passed"},
            "sources": [record(csrc / (name + ".c")) for name in KERNELS] + [record(root / "setup.py")],
            "binaries": [record(csrc / (name + ".so")) for name in KERNELS]}
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(data))
    return root, manifest, data


def test_runtime_complete(tmp_path):
    root, manifest, _ = fixture_manifest(tmp_path)
    assert verify_runtime(manifest, root)["status"] == "complete"


@pytest.mark.parametrize("suffix", [".so", ".c"])
def test_changed_file_rejected(tmp_path, suffix):
    root, manifest, _ = fixture_manifest(tmp_path)
    (root / ("orthohmm/search/csrc/pair_align" + suffix)).write_text("changed")
    with pytest.raises(ValueError, match="Changed native"):
        verify_runtime(manifest, root)


def test_extra_cuda_rejected(tmp_path):
    root, manifest, _ = fixture_manifest(tmp_path)
    (root / "orthohmm/search/csrc/hmm_viterbi_cuda.so").write_bytes(b"cuda")
    with pytest.raises(ValueError, match="library set"):
        verify_runtime(manifest, root)


@pytest.mark.parametrize("field,value", [("status", "failed"), ("commit", "wrong"),
                                        ("profile_probe", {}), ("binaries", []), ("sources", [])])
def test_incomplete_manifest_rejected(tmp_path, field, value):
    root, manifest, data = fixture_manifest(tmp_path)
    data[field] = value
    manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        verify_runtime(manifest, root)
