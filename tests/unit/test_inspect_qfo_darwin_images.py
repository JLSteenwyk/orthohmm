import io
import tarfile

import pytest

from benchmark_tools.inspect_qfo_darwin_images import context_layers, digest, inventory, verified


def test_digest_and_size_both_required():
    data = b"public blob"
    descriptor = {"size": len(data), "digest": digest(data)}
    assert verified(data, descriptor) == data
    with pytest.raises(ValueError):
        verified(data, {**descriptor, "size": len(data) + 1})
    with pytest.raises(ValueError):
        verified(data + b"!", {**descriptor, "size": len(data) + 1})


def test_history_skips_empty_entries_and_maps_copy_only():
    layers = [{"digest": "base"}, {"digest": "context"}, {"digest": "other"}]
    history = [{"created_by": "ADD ubuntu.tar /"},
               {"created_by": "WORKDIR /benchmark", "empty_layer": True},
               {"created_by": "COPY . /benchmark # buildkit"},
               {"created_by": "RUN install"}]
    assert context_layers({"layers": layers}, {"history": history}) == [
        (1, layers[1], "COPY . /benchmark # buildkit")]


def test_bad_history_does_not_imply_no_originals():
    with pytest.raises(ValueError):
        context_layers({"layers": [{}]}, {"history": []})


def test_tar_is_inventory_only_even_with_unsafe_paths_and_links():
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name in ["benchmark/data/treefam/TF0001.nhx.gz",
                     "../treefam2reference.txt", "benchmark/TreeFam-A.json"]:
            member = tarfile.TarInfo(name)
            member.size = 3
            archive.addfile(member, io.BytesIO(b"abc"))
        link = tarfile.TarInfo("benchmark/external")
        link.type = tarfile.SYMTYPE
        link.linkname = "/etc/passwd"
        archive.addfile(link)
    result = inventory(output.getvalue())
    assert result["member_count"] == 4
    assert len(result["original_filename_candidates"]) == 2
    assert len(result["treefam_named_members"]) == 3
    assert result["members"][-1]["linkname"] == "/etc/passwd"
