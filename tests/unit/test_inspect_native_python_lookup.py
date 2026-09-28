import pytest
from copy import deepcopy

from benchmark_tools.inspect_native_python_lookup import covered_files, scientific_origins, compare_lookup


def test_direct_and_symlink_file_coverage():
    rows = [dict(kind="file", path="/a", bytes=4, sha256="a"),
            dict(kind="symlink", path="/b", resolved="/target", target_bytes=5, target_sha256="b")]
    files = [dict(path="/a", bytes=4, sha256="a"), dict(path="/target", bytes=5, sha256="b")]
    assert covered_files(files, rows) == dict(missing=[], changed=[], all_covered=True)
    files[0]["sha256"] = "changed"
    files.append(dict(path="/new", bytes=0, sha256="new"))
    assert covered_files(files, rows) == dict(missing=["/new"], changed=["/a"], all_covered=False)


def test_scientific_source_root_not_merely_inventory_coverage():
    report = dict(modules={"orthohmm": "/frozen/orthohmm/__init__.py",
                           "orthohmm.search": "/frozen/orthohmm/search/__init__.py"})
    assert scientific_origins(report, "orthohmm", "/frozen/orthohmm")["checked_modules"] == 2
    report["modules"]["orthohmm.search"] = "/frozen/orthohmm_other/search.py"
    with pytest.raises(ValueError):
        scientific_origins(report, "orthohmm", "/frozen/orthohmm")
    with pytest.raises(ValueError):
        scientific_origins(report, "orthofinder", "/frozen/orthofinder")


@pytest.mark.parametrize("field", ["paths", "modules", "mapped_files", "meta_path", "path_hooks", "editable", "files"])
def test_lookup_drift_is_rejected(field):
    expected = dict(python="version", executable="/python", cwd="/core", requested=["core"],
        paths=[], modules={}, mapped_files=[], meta_path=[], path_hooks=[], editable={}, files=[],
        dont_write_bytecode=True, coverage=dict(missing=[], changed=[], all_covered=True),
        scientific_origin=dict(root="/core"), pycache_prefix="/fresh/a", ipc_mappings=[])
    observed = deepcopy(expected)
    observed.update(pycache_prefix="/fresh/b", ipc_mappings=[dict(path="/dev/shm/sem.example")])
    assert compare_lookup(expected, observed)["status"] == "native_python_lookup_matches"
    observed[field] = ["changed"]
    with pytest.raises(ValueError):
        compare_lookup(expected, observed)
