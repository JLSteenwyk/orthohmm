import json

import pytest

import benchmark_tools.audit_cpm_memcheck_artifacts as module


def setup(tmp_path, monkeypatch):
    directory = tmp_path / "benchmarks/results/qfo_cpm_refinement_memcheck_diagnostic_v1"
    native = tmp_path / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    directory.mkdir(parents=True)
    (native / "payload").mkdir(parents=True)
    names = native / "payload/gene_names.txt"
    names.write_text("a\nb\n")
    output = directory / "refinement_repeat.txt"
    output.write_text("a b\n")
    native_output = native / "groups.txt"
    native_output.write_text("a b\n")
    metadata = dict(genes=2, groups=1, accuracy_evaluated=False)
    (native / "refinement.json").write_text(json.dumps(dict(metadata, output=module.record(native_output))))
    (directory / "refinement_repeat.json").write_text(json.dumps(dict(metadata, output=module.record(output))))
    log, xml = directory / "child.log", directory / "memcheck.xml"
    log.write_text("diagnostic failed")
    xml.write_text("retained XML placeholder; parsing tested separately")
    failed = dict(status="memcheck_diagnostic_failed", returncode=97,
                  checked_records=[module.record(names), module.record(native / "refinement.json")],
                  runtime_before={}, child_log=module.record(log), memcheck_xml=module.record(xml))
    status = directory / "status.json"
    status.write_text(json.dumps(failed))
    monkeypatch.setattr(module, "STATUS_SHA", module.record(status)["sha256"])
    return directory, native


def test_checks_saved_partition_without_admitting_failure(tmp_path, monkeypatch):
    setup(tmp_path, monkeypatch)
    result = module.audit(tmp_path)
    assert result["coverage"] == dict(genes=2, groups=1)
    assert result["original_returncode"] == 97
    assert result["before"]["mismatches"] == result["after"]["mismatches"] == []
    assert not result["seed_admitted"] and not result["publication_ready"]


@pytest.mark.parametrize("target", ["child.log", "refinement_repeat.txt", "refinement_repeat.json"])
def test_changed_artifact_rejected(tmp_path, monkeypatch, target):
    directory, native = setup(tmp_path, monkeypatch)
    path = directory / target
    if target.endswith(".json"):
        data = json.loads(path.read_text())
        data["groups"] = 9
        path.write_text(json.dumps(data))
    else:
        path.write_text("changed")
    with pytest.raises(ValueError):
        module.audit(tmp_path)


def test_identity_inventory_rejects_conflicts_and_records_missing(tmp_path):
    path = tmp_path / "input"
    path.write_text("same")
    item = module.record(path)
    assert module.inspect([item, item])["unique_files"] == 1
    with pytest.raises(ValueError, match="Conflicting"):
        module.inspect([item, dict(item, bytes=3)])
    path.unlink()
    assert module.inspect([item])["matching_files"] == 0


def test_collect_nested_runtime_records():
    item = dict(path="x", bytes=1, sha256="a")
    assert list(module.collect(dict(runtime=[dict(files=[item])], count=1))) == [item]
