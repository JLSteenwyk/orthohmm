import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

from tests.unit import test_bundle_publication_handoff as legacy
from tests.unit import test_bundle_publication_source as source_legacy

candidate = legacy.candidate
native_candidate = legacy.native_candidate
exported = source_legacy.exported
MAIN = "benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261004.md"
STAGES = {role: "benchmark_tools/results/final_" + role + ".json"
          for role in ("render", "print", "review")}


def chosen(revision="a" * 40):
    return dict(review_revision=revision, ledger_revision=revision,
                main_text=MAIN, stages=dict(STAGES))


def refresh(directory, index):
    module = legacy.module
    for role, (index_name, _) in module.COMPONENTS.items():
        path = directory / role / index_name
        index["components"][role]["sha256"] = module.identity(path.read_bytes())["sha256"]
        index["components"][role]["verified"]["manifest"] = module.identity(path.read_bytes())
    index["files"] = [dict(path=p.relative_to(directory).as_posix(),
        mode=p.stat().st_mode & 0o777, **module.identity(p.read_bytes()))
        for p in sorted(directory.rglob("*")) if p.is_file()
        and p != directory / "HANDOFF_INDEX.json"]
    path = directory / "HANDOFF_INDEX.json"
    path.write_text(json.dumps(index))
    return module.identity(path.read_bytes())["sha256"]


@pytest.fixture(params=["native-preparation", "native-build"])
def selected(request):
    module = legacy.module
    directory, _ = request.getfixturevalue("candidate" if request.param == "native-preparation" else "native_candidate")
    index = json.loads((directory / "HANDOFF_INDEX.json").read_bytes())
    index.update(schema="publication_handoff_candidate_v3", source_profile=request.param,
                 result_helpers_included=True, review_selection=dict(chosen(), page_count=12))
    for role, (index_name, verifier_name) in module.COMPONENTS.items():
        path = directory / role / index_name
        child = json.loads(path.read_bytes())
        child["files"][0]["git_revision"] = "a" * 40
        result = dict(publication_ready=False, page_count=12, fixture_interface=True)
        if role == "source":
            child.update(schema="publication_source_components_v3", workflow_revision="a" * 40,
                         include_result_helpers=True)
            result["result_helpers_included"] = True
        else:
            child.update(schema="publication_direct_review_v3", main_text=MAIN, stages=dict(STAGES))
            result["entrypoints"] = dict(markdown=MAIN)
        path.write_text(json.dumps(child))
        result["manifest"] = module.identity(path.read_bytes())
        index["components"][role]["verified"] = result
        verifier = directory / role / verifier_name
        # Synthetic interfaces exercise coordinator binding, not actual child exports.
        code = ("import hashlib\nfrom pathlib import Path\n"
            f"RUNNER={verifier_name!r}\nGUIDE='fixture-guide.md'\n"
            "LEDGER='benchmark_tools/results/PUBLICATION_PROGRESS.md'\n"
            "def manuscript_path(name):\n    return name\n"
            "def stage_paths(value,main):\n    return value\n"
            "def verify(directory,digest):\n"
            f"    content=(Path(directory)/{index_name!r}).read_bytes()\n"
            "    if hashlib.sha256(content).hexdigest()!=digest:\n        raise ValueError('Component digest differs')\n"
            f"    result={result!r}\n"
            "    result['manifest']=dict(bytes=len(content),sha256=digest)\n"
            "    return result\n")
        verifier.write_text(code)
    return directory, refresh(directory, index)


def test_explicit_review_and_helpers_relocated(selected, tmp_path):
    module = legacy.module
    directory, digest = selected
    result = module.verify(directory, digest)
    assert result["review_selection"] == dict(chosen(), page_count=12)
    assert result["result_helpers_included"] is True
    assert result["components"]["manuscript"]["page_count"] == 12
    assert result["publication_ready"] is False
    copied = tmp_path / "relocated-explicit-review"
    shutil.move(directory, copied)
    checked = subprocess.run([sys.executable, "-I", "-S", "-B",
        str(copied / "bundle_publication_handoff.py"), "verify", str(copied),
        "--manifest-sha256", digest], cwd=tmp_path,
        env={"PATH": "/no-git-or-original-artifacts"}, capture_output=True,
        text=True, timeout=30)
    assert checked.returncode == 0, checked.stderr
    assert json.loads(checked.stdout) == result


@pytest.mark.parametrize("defect", ["missing_selection", "helper_flag", "legacy_v1", "legacy_v2",
    "page_bool", "page_mismatch", "main_mismatch", "stage_mismatch", "source_helpers",
    "source_revision", "review_revision", "review_schema"])
def test_explicit_component_bindings_not_relaxed(selected, defect):
    module = legacy.module
    directory, _ = selected
    index = json.loads((directory / "HANDOFF_INDEX.json").read_bytes())
    if defect == "missing_selection":
        del index["review_selection"]
    elif defect == "helper_flag":
        index["result_helpers_included"] = False
    elif defect.startswith("legacy_"):
        index["schema"] = "publication_handoff_candidate_" + defect.split("_", 1)[1]
        if defect == "legacy_v1":
            index.pop("source_profile")
    elif defect.startswith("page_"):
        index["review_selection"]["page_count"] = True if defect == "page_bool" else 13
    else:
        role = "source" if defect.startswith("source_") else "manuscript"
        path = directory / role / module.COMPONENTS[role][0]
        child = json.loads(path.read_bytes())
        if defect == "main_mismatch": child["main_text"] = "wrong.md"
        elif defect == "stage_mismatch": child["stages"]["render"] = "wrong.json"
        elif defect == "source_helpers": child["include_result_helpers"] = False
        elif defect == "source_revision": child["workflow_revision"] = "b" * 40
        elif defect == "review_revision": child["files"][0]["git_revision"] = "b" * 40
        else: child["schema"] = "publication_direct_review_v2"
        path.write_text(json.dumps(child))
    digest = refresh(directory, index)
    with pytest.raises(ValueError):
        module.verify(directory, digest)


@pytest.mark.parametrize("change", [
    {"main_text": "../main.md"}, {"main_text": "main.json"},
    {"review_revision": "HEAD"}, {"ledger_revision": True}, {"extra": 1},
    {"stages": {}}, {"stages": dict(render=[], print="b.json", review="c.json")},
    {"stages": dict(render="a.json", print="a.json", review="c.json")},
])
def test_bad_selection_fails_before_output(tmp_path, change):
    selection = chosen()
    selection.update(change)
    with pytest.raises(ValueError):
        legacy.module.build(tmp_path, "HEAD", tmp_path / "output", selected_review=selection)
    assert not (tmp_path / "output").exists()


def test_builder_forwards_explicit_selection_and_records_it(tmp_path, monkeypatch):
    module = legacy.module
    import benchmark_tools
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    for key, value in (("user.name", "Test"), ("user.email", "test@example.invalid")):
        subprocess.run(["git", "-C", str(repo), "config", key, value], check=True)
    fixture_root = Path(module.__file__).resolve().parent.parent
    for name, _ in module.EXTRAS.values():
        target = repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((fixture_root / name).read_bytes())
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "handoff fixture"], check=True)
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    calls = []

    def source_build(repo, commit, output, profile, **options):
        calls.append(("source", commit, profile, options))
        output.mkdir()
        (output / "SOURCE_INDEX.json").write_text("{}")
        return dict(manifest=module.identity(b"{}"), page_count=12)

    def review_build(repo, review, ledger, workflow, output, **options):
        calls.append(("review", review, ledger, workflow, options))
        output.mkdir()
        (output / "REVIEW_INDEX.json").write_text("{}")
        return dict(manifest=module.identity(b"{}"), page_count=12)

    monkeypatch.setattr(benchmark_tools, "bundle_publication_source", SimpleNamespace(build=source_build))
    committed = benchmark_tools.bundle_publication_review.committed
    monkeypatch.setattr(benchmark_tools, "bundle_publication_review", SimpleNamespace(build=review_build,
        committed=committed, manuscript_path=lambda main: main, stage_paths=lambda stages, main: stages))
    monkeypatch.setattr(module, "verify", lambda directory, digest:
                        json.loads((directory / "HANDOFF_INDEX.json").read_bytes()))
    result = module.build(repo, revision, tmp_path / "export", selected_review=chosen(revision))
    assert calls == [("source", revision, "native-preparation", dict(include_result_helpers=True)),
                     ("review", revision, revision, revision, dict(main_text=MAIN, stages=STAGES))]
    assert result["schema"] == "publication_handoff_candidate_v3"
    assert result["review_selection"] == dict(chosen(revision), page_count=12)
    assert result["result_helpers_included"] is True
    assert "nine-page" not in " ".join(result["limitations"])


def test_cli_requires_complete_explicit_review(tmp_path):
    completed = subprocess.run([sys.executable, "-I", "-S", "-B", legacy.module.__file__,
        "build", "--repo", str(tmp_path), "--revision", "HEAD", "--output", str(tmp_path / "output"),
        "--main-text", MAIN], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert completed.returncode == 2
    assert "requires both revisions" in completed.stderr
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("exported", ["native-preparation", "native-build"], indirect=True)
def test_actual_child_builders_and_relocated_handoff(exported, tmp_path):
    import benchmark_tools
    module = legacy.module
    review = benchmark_tools.bundle_publication_review
    repo, _, _ = exported
    fixture_root = Path(module.__file__).resolve().parent.parent
    profile = json.loads((exported[1] / "SOURCE_INDEX.json").read_bytes())["profile"]
    for name, _ in module.selection(profile).values():
        target = repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((fixture_root / name).read_bytes())
    html = "benchmark_tools/results/final_review.html"
    pdf = "benchmark_tools/results/final_print/document.pdf"
    page = "benchmark_tools/results/final_pdf/page_001.png"
    payloads = {MAIN: b"# Final fixture\n[Ledger](PUBLICATION_PROGRESS.md)\n",
        review.LEDGER: b"Frozen fixture ledger\n",
        html: b'<html><a href="PUBLICATION_PROGRESS.md">Ledger</a></html>',
        pdf: b"Synthetic PDF fixture", page: b"Synthetic page fixture",
        review.RUNNER: Path(review.__file__).read_bytes(), review.GUIDE: b"Review fixture scope",
        "benchmark_tools/results/fixture_helper.py": b"def marker():\n    return 42\n"}
    for name, content in payloads.items():
        target = repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)

    def record(name):
        return dict(path="/historical/worktree/" + name,
                    **review.identity((repo / name).read_bytes()))

    def save(name, data):
        target = repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(data))

    save(STAGES["render"], dict(status="manuscript_review_rendered", publication_ready=False,
        html=record(html), sources=[record(MAIN)], targets=[record(review.LEDGER)], unique_targets=1,
        occurrences=[dict(path=review.LEDGER, url="PUBLICATION_PROGRESS.md")], local_occurrences=1))
    save(STAGES["print"], dict(status="verified_html_printed", publication_ready=False,
        returncode=0, page_count=1, pdf=record(pdf), checked_records=[record(STAGES["render"]), record(html)]))
    save(STAGES["review"], dict(status="pdf_bounds_checked", publication_ready=False,
        page_count=1, bounds_violations=[], checked_records=[record(pdf), record(STAGES["render"])],
        rendered_pages=[record(page)]))
    # The source fixture intentionally dirties probe.py; do not stage that file.
    paths = {*[name for name, _ in module.selection(profile).values()], *payloads, *STAGES.values()}
    subprocess.run(["git", "-C", str(repo), "add", *sorted(paths)], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "actual child fixture integration"], check=True)
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    output = tmp_path / "combined-export"
    result = module.build(repo, revision, output, profile, selected_review=chosen(revision))
    assert result["result_helpers_included"] is True
    assert result["review_selection"]["page_count"] == 1
    assert result["components"]["manuscript"]["entrypoints"]["markdown"] == MAIN
    assert result["publication_ready"] is False
    assert (output / "source/workflow/benchmark_tools/results/fixture_helper.py").exists()
    assert (output / "source/scientific/orthohmm/version.py").read_text() == "VERSION = 'fixture'\n"
    relocated = tmp_path / "relocated-real-children"
    shutil.move(output, relocated)
    shutil.rmtree(repo)
    checked = subprocess.run([sys.executable, "-I", "-S", "-B",
        str(relocated / "bundle_publication_handoff.py"), "verify", str(relocated),
        "--manifest-sha256", result["manifest"]["sha256"]], cwd=tmp_path,
        env={"PATH": "/no-git-or-original-artifacts"}, capture_output=True,
        text=True, timeout=30)
    assert checked.returncode == 0, checked.stderr
    assert json.loads(checked.stdout) == result
