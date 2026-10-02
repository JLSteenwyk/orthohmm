import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import bundle_publication_handoff as module


@pytest.fixture
def candidate(tmp_path):
    root = Path(module.__file__).resolve().parent.parent
    directory = tmp_path / "handoff"
    directory.mkdir()
    extras = {}
    for name, (git_path, _) in module.EXTRAS.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((root / git_path).read_bytes())
        path.chmod(0o644)
        extras[name] = dict(git_path=git_path, git_revision="a" * 40, git_blob="b" * 40)
    components = {}
    # Synthetic component interfaces isolate coordinator checks; real component
    # verifiers and the actual committed export are exercised separately.
    for role, (index_name, verifier_name) in module.COMPONENTS.items():
        component = directory / role
        verifier = component / verifier_name
        verifier.parent.mkdir(parents=True)
        code = ("import hashlib\nfrom pathlib import Path\n"
                "def verify(directory, digest):\n"
                f"    content = (Path(directory) / {index_name!r}).read_bytes()\n"
                "    if hashlib.sha256(content).hexdigest() != digest:\n"
                "        raise ValueError('Component digest differs')\n"
                "    return dict(publication_ready=False, page_count=9, fixture_interface=True)\n")
        verifier.write_text(code)
        verifier.chmod(0o644)
        child = dict(profile="native-preparation" if role == "source" else None,
            files=[dict(path=verifier_name)])
        index = component / index_name
        index.write_text(json.dumps(child))
        index.chmod(0o644)
        components[role] = dict(index=role + "/" + index_name,
            sha256=module.identity(index.read_bytes())["sha256"],
            verified=dict(publication_ready=False, page_count=9, fixture_interface=True))
    files = [dict(path=p.relative_to(directory).as_posix(), mode=0o644, **module.identity(p.read_bytes()))
             for p in sorted(directory.rglob("*")) if p.is_file()]
    index = dict(schema="publication_handoff_candidate_v1", files=files, components=components,
        extra_sources=extras, workflow_revision="a" * 40, publication_ready=False,
        redistribution_clearance=False, public_release_uploaded=False)
    (directory / "HANDOFF_INDEX.json").write_text(json.dumps(index))
    return directory, module.identity((directory / "HANDOFF_INDEX.json").read_bytes())["sha256"]


def update_index(directory, change):
    path = directory / "HANDOFF_INDEX.json"
    data = json.loads(path.read_bytes())
    change(data)
    path.write_text(json.dumps(data))
    return module.identity(path.read_bytes())["sha256"]


def test_candidate_components_and_scope(candidate):
    directory, digest = candidate
    result = module.verify(directory, digest)
    assert result["comparison"] == dict(methods=8, numeric_cells=72,
        raw_scoring_repeated=False, secondary_mean_official=False)
    assert set(result["components"]) == {"source", "manuscript"}
    assert result["publication_ready"] is False
    assert result["native_inference_reproduced"] is False
    assert result["numerical_replay_executed"] is False


def test_relocated_cli_without_original_directory_or_git(candidate, tmp_path):
    directory, digest = candidate
    copied = tmp_path / "copied"
    shutil.copytree(directory, copied)
    shutil.rmtree(directory)
    result = subprocess.run([sys.executable, "-I", "-S", "-B",
        str(copied / "bundle_publication_handoff.py"), "verify", str(copied),
        "--manifest-sha256", digest], cwd=tmp_path, env={"PATH": "/no-git-here"},
        capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["comparison"]["numeric_cells"] == 72


@pytest.mark.parametrize("name", ["", ".", "../a", "/a", "a//b", "a/./b", "a/../b", "a\\b", 1])
def test_unsafe_paths(name):
    with pytest.raises(ValueError):
        module.relative(name)


@pytest.mark.parametrize("defect", ["payload", "missing", "extra", "symlink", "mode", "digest"])
def test_reject_inventory_before_component_code(candidate, monkeypatch, defect):
    directory, digest = candidate
    path = directory / "README.md"
    if defect == "payload":
        path.write_text("changed")
    elif defect == "missing":
        path.unlink()
    elif defect == "extra":
        (directory / "extra").write_text("extra")
    elif defect == "symlink":
        path.unlink()
        path.symlink_to(directory / "LICENSE.md")
    elif defect == "mode":
        path.chmod(0o777)
    else:
        digest = "0" * 64
    monkeypatch.setattr(module, "load", lambda *args: pytest.fail("Do not load before inventory verification"))
    with pytest.raises((ValueError, FileNotFoundError)):
        module.verify(directory, digest)


@pytest.mark.parametrize("defect", ["duplicate", "ready", "clearance", "upload", "profile", "mapping", "component_result"])
def test_changed_index_rejected_with_reanchored_digest(candidate, defect):
    directory, _ = candidate
    def change(index):
        if defect == "duplicate":
            index["files"].append(index["files"][0])
        elif defect in {"ready", "clearance", "upload"}:
            index[{"ready": "publication_ready", "clearance": "redistribution_clearance",
                   "upload": "public_release_uploaded"}[defect]] = True
        elif defect == "profile":
            index["components"]["source"]["index"] = "../SOURCE_INDEX.json"
        elif defect == "mapping":
            index["extra_sources"]["README.md"]["git_revision"] = "f" * 40
        else:
            index["components"]["source"]["verified"]["page_count"] = 8
    digest = update_index(directory, change)
    with pytest.raises(ValueError):
        module.verify(directory, digest)


def test_coordinated_content_index_change_requires_external_anchor(candidate):
    directory, digest = candidate
    path = directory / "README.md"
    path.write_text("changed")
    def change(index):
        for row in index["files"]:
            if row["path"] == "README.md":
                row.update(module.identity(path.read_bytes()))
    update_index(directory, change)
    with pytest.raises(ValueError, match="external anchor"):
        module.verify(directory, digest)


def test_existing_build_output_refused(candidate):
    directory, _ = candidate
    with pytest.raises(FileExistsError):
        module.build(directory, "HEAD", directory)


@pytest.mark.parametrize("value", ["nan", "inf", "-1", "2"])
def test_table_invalid_numeric_values(candidate, value):
    directory, _ = candidate
    content = (directory / "comparison/scores.tsv").read_text()
    content = content.replace("0.7035899765495244", value)
    with pytest.raises(ValueError):
        module.check_table(content.encode())


def test_table_secondary_mean_is_not_recomputed_as_f1(candidate):
    directory, _ = candidate
    content = (directory / "comparison/scores.tsv").read_text().replace("0.6895554410913292", "0.9")
    with pytest.raises(ValueError, match="secondary mean differs"):
        module.check_table(content.encode())
