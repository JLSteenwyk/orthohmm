from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import pytest

from benchmark_tools import bundle_qfo_parameter_component as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_bootstrap_qfo_parameter_neighborhood import fixture as counts_fixture
from tests.unit.test_reproduce_qfo_parameter_uncertainty import result_for


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")
    return record(path)


def commit(repo):
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
        "-c", "commit.gpgsign=false", "commit", "-qm", "fixture"], check=True)


@pytest.fixture
def repository(tmp_path):
    root = Path(module.__file__).resolve().parents[1]
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    for name in (module.RUNNER, "LICENSE.md", *["benchmark_tools/" + name for name in module.SOURCES]):
        target = repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / name, target)
    analysis = result_for(counts_fixture())
    analysis.update(controlled_timing=False, source=record(repo / "benchmark_tools/run_private_cpm_parameter_uncertainty.py"),
        helpers=[record(repo / "benchmark_tools" / name) for name in module.SOURCES if name != "run_private_cpm_parameter_uncertainty.py"])
    for key, name in (("protocol", "scientific.md"), ("execution_protocol", "execution.md"), ("plan", "plan.json")):
        target = repo / "benchmark_tools/results" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("fixture metadata; not a scientific protocol\n")
        analysis[key] = record(target)
    pin = save(repo / "benchmarks/work/qfo_private_cpm_parameter_uncertainty_20261001.json", analysis)
    reproduction = dict(status="qfo_parameter_uncertainty_numerically_reproduced", input=pin, endpoints=18,
        planned_endpoints=18, absolute_tolerance=1e-12, publication_ready=False,
        numpy_version=analysis["numpy_version"], source=record(repo / "benchmark_tools" / module.REPRODUCER))
    repro_pin = save(repo / "benchmark_tools/results/qfo_private_cpm_parameter_reproduction_20261001.json", reproduction)
    summary = dict(status="qfo_private_cpm_parameter_result_summary", source=analysis["source"], input=pin,
        reproduction=repro_pin, point_estimates=analysis["point_estimates"], families=analysis["families"],
        comparisons=[{key: row[key] for key in ("candidate", "reference", "status", "metrics")} for row in analysis["comparisons"]],
        complete_panel=True, endpoint_count=18, estimated_contrasts=6, replicates=100000, seed=20260925,
        multiplicity_endpoints=18, publication_ready=False, controlled_timing=False)
    save(repo / module.SUMMARY, summary)
    outputs = []
    for name in module.FIGURE_NAMES:
        target = repo / "benchmark_tools/results/qfo_parameter_complete_export_20261001" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("test output placeholder; not rendered media or a real metric table\n")
        outputs.append(record(target))
    exported = dict(status="qfo_parameter_results_exported", input=pin, reproduction=repro_pin,
        outputs=outputs, complete_panel=True, historical_exporter_binding=None, estimated_endpoints=18,
        planned_endpoints=18, publication_ready=False)
    export_pin = save(repo / "benchmark_tools/results/qfo_parameter_complete_export_20261001/manifest.json", exported)
    save(repo / module.READBACK, dict(status="qfo_parameter_export_output_readback", input=pin,
        reproduction=repro_pin, outputs=outputs, complete_panel=True, historical_exporter_binding=None,
        estimated_endpoints=18, planned_endpoints=18, publication_ready=False, manifest=export_pin))
    commit(repo)
    return repo


def built(repository, tmp_path):
    output = tmp_path / "component"
    result = module.build(repository, "HEAD", output)
    return output, result["manifest"]["sha256"], result


def test_relocation_and_actual_isolated_numerical_cli_without_original_checkout(repository, tmp_path):
    output, sha, expected = built(repository, tmp_path)
    moved = tmp_path / "relocated"
    output.rename(moved)
    repository.rename(tmp_path / "unavailable-original-checkout")
    assert module.verify(moved, sha)[0] == expected
    run = subprocess.run([sys.executable, "-I", "-B", str(moved / module.RUNNER), "reproduce",
        "--directory", str(moved), "--manifest-sha256", sha, "--output", str(tmp_path / "numeric.json")],
        cwd=tmp_path, capture_output=True, text=True, check=True)
    receipt = json.loads(run.stdout)
    assert receipt["status"] == "qfo_parameter_component_numerically_reproduced"
    assert receipt["endpoints"] == 18 and receipt["absolute_tolerance"] == 1e-12
    assert receipt["native_inference"] is False and receipt["scientific_inputs_readmitted"] is False
    assert receipt["publication_ready"] is False
    assert json.loads((tmp_path / "numeric.json").read_text()) == receipt


@pytest.mark.parametrize("damage", ["bytes", "missing", "extra", "symlink", "directory_symlink", "manifest_symlink",
    "manifest", "traversal", "duplicate", "scope", "claims", "version"])
def test_damaged_or_promoted_component_rejected(repository, tmp_path, damage):
    output, sha, _ = built(repository, tmp_path)
    path = output / "figures/scores.tsv"
    manifest_path = output / "bundle.json"
    manifest = json.loads(manifest_path.read_text())
    if damage == "bytes":
        path.write_text("changed")
    elif damage == "missing":
        path.unlink()
    elif damage == "extra":
        (output / "unexpected").write_text("extra")
    elif damage == "symlink":
        target = tmp_path / "source.tsv"
        shutil.copyfile(path, target)
        path.unlink()
        path.symlink_to(target)
    elif damage == "directory_symlink":
        (output / "extra_directory").symlink_to(repository, target_is_directory=True)
    elif damage == "manifest_symlink":
        target = tmp_path / "manifest.json"
        manifest_path.rename(target)
        manifest_path.symlink_to(target)
    elif damage == "manifest":
        manifest["numpy_version"] = "wrong"
        manifest_path.write_text(json.dumps(manifest))
    else:
        if damage == "traversal":
            manifest["files"][0]["path"] = "../escape"
        elif damage == "duplicate":
            manifest["files"].append(manifest["files"][0])
        elif damage == "scope":
            manifest["scope"] = "full_native_reproduction"
        elif damage == "claims":
            manifest["publication_ready"] = True
        else:
            manifest["numpy_version"] = "wrong"
        manifest_path.write_text(json.dumps(manifest))
        sha = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    with pytest.raises((ValueError, FileNotFoundError)):
        module.verify(output, sha)


def test_committed_support_not_dirty_checkout_and_changed_local_only_evidence_rejected(repository, tmp_path):
    (repository / "LICENSE.md").write_text("dirty local license bytes")
    output, sha, _ = built(repository, tmp_path)
    assert (output / "LICENSE.md").read_text() != "dirty local license bytes"
    (repository / "benchmarks/work/qfo_private_cpm_parameter_uncertainty_20261001.json").write_text("changed")
    with pytest.raises(ValueError, match="evidence bytes"):
        module.build(repository, "HEAD", tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_failed_reproduction_receipt_retained_and_no_overwrite(repository, tmp_path):
    output, sha, _ = built(repository, tmp_path)
    (output / "data/analysis.json").write_text("changed")
    receipt = tmp_path / "failure.json"
    with pytest.raises(ValueError):
        module.reproduce(output, sha, receipt)
    assert json.loads(receipt.read_text())["status"] == "validation_failed"
    with pytest.raises(FileExistsError):
        module.reproduce(output, sha, receipt)


def test_executing_builder_must_match_selected_git_revision(repository, tmp_path):
    with (repository / module.RUNNER).open("a") as stream:
        stream.write("\n# changed committed builder\n")
    commit(repository)
    with pytest.raises(ValueError, match="Executing builder differs"):
        module.build(repository, "HEAD", tmp_path / "wrong_source")
    assert not (tmp_path / "wrong_source").exists()


def test_declared_numpy_version_is_required_before_arithmetic(repository, tmp_path, monkeypatch):
    output, sha, _ = built(repository, tmp_path)
    function, np = module.numerical_verifier((output / "sources" / module.REPRODUCER).read_bytes())
    monkeypatch.setattr(np, "__version__", "not-the-declared-version")
    with pytest.raises(ValueError, match="declared numerical"):
        module.reproduce(output, sha, tmp_path / "wrong_numpy.json")
    assert json.loads((tmp_path / "wrong_numpy.json").read_text())["status"] == "validation_failed"


@pytest.mark.parametrize("damage", ["summary_point", "summary_complete", "summary_source", "source_digest",
    "reproduction_endpoint", "reproduction_input", "export_scope", "export_output_order", "protocol_bytes"])
def test_semantic_binding_rejects_coherently_rewrapped_metadata(repository, tmp_path, damage):
    output, sha, _ = built(repository, tmp_path)
    _, payloads = module.verify(output, sha)
    if damage == "protocol_bytes":
        payloads["data/scientific_protocol.md"] = b"changed protocol"
    else:
        name = "data/summary.json" if damage.startswith("summary") else (
            "data/analysis.json" if damage == "source_digest" else
            "data/reproduction.json" if damage.startswith("reproduction") else "data/export_manifest.json")
        data = json.loads(payloads[name])
        if damage == "summary_point":
            data["point_estimates"]["control"]["F1"] += .1
        elif damage == "summary_complete":
            data["complete_panel"] = False
        elif damage in ("summary_source", "source_digest"):
            data["source"]["sha256"] = "changed"
        elif damage == "reproduction_endpoint":
            data["endpoints"] = 15
        elif damage == "reproduction_input":
            data["input"]["sha256"] = "changed"
        elif damage == "export_scope":
            data["publication_ready"] = True
        else:
            data["outputs"].reverse()
        payloads[name] = json.dumps(data).encode()
    with pytest.raises(ValueError):
        module.semantic_checks(payloads)


def test_archive_is_deterministic_all_regular_and_restores_from_fresh_extraction(repository, tmp_path):
    output, sha, _ = built(repository, tmp_path)
    a, b = tmp_path / "first.tar.gz", tmp_path / "second.tar.gz"
    result = module.archive(output, sha, a)
    assert module.archive(output, sha, b)["archive"] == result["archive"]
    assert a.read_bytes() == b.read_bytes()
    restored = tmp_path / "extracted"
    restored.mkdir()
    with tarfile.open(a, "r:gz") as handle:
        for member in handle.getmembers():
            assert member.isfile()
            path = restored / module.safe_name(member.name)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(handle.extractfile(member).read())
    assert module.verify(restored, sha)[0]["manifest"]["sha256"] == sha


def test_fresh_and_outside_component_outputs_required(repository, tmp_path):
    output, sha, _ = built(repository, tmp_path)
    with pytest.raises(FileExistsError):
        module.build(repository, "HEAD", output)
    with pytest.raises(ValueError):
        module.reproduce(output, sha, output / "result.json")
    with pytest.raises(ValueError):
        module.archive(output, sha, output / "bundle.tar.gz")
    target = tmp_path / "existing.tar.gz"
    target.write_bytes(b"retained")
    with pytest.raises(FileExistsError):
        module.archive(output, sha, target)
    assert target.read_bytes() == b"retained"


@pytest.mark.parametrize("damage", ["counts", "point", "interval", "truth", "order", "source"])
def test_original_pure_verifier_is_reused_and_rejects_numerical_corruption(damage):
    root = Path(module.__file__).resolve().parents[1]
    source = (root / "benchmark_tools" / module.REPRODUCER).read_bytes()
    if damage == "source":
        with pytest.raises(ValueError):
            module.numerical_verifier(source + b"\n# changed\n")
        return
    report = deepcopy(result_for(counts_fixture()))
    if damage == "counts":
        report["reconstructed_counts"]["arms"][0]["families"][0]["counts_without_prior"]["TP"] = -1
    elif damage == "truth":
        report["reconstructed_counts"]["arms"][1]["families"][0]["counts_without_prior"]["TN"] += 1
    elif damage == "order":
        report["comparisons"].reverse()
    else:
        item = report["comparisons"][0]["metrics"]["F1"]
        if damage == "point":
            item["difference"] += .01
        else:
            item["paired_percentile_ci"][0] += .01
    function, _ = module.numerical_verifier(source)
    with pytest.raises(ValueError):
        function(report)
