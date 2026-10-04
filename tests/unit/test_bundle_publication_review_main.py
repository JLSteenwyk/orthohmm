import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from tests.unit import test_bundle_publication_review as legacy

prepared = legacy.prepared
MAIN = "benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261004.md"


@pytest.fixture
def selected(prepared):
    module = legacy.bundle
    repo, _, _, output = prepared
    old = repo / module.MAIN
    new = repo / MAIN
    new.write_bytes(old.read_bytes())
    old.unlink()
    render = json.loads((repo / module.RENDER).read_bytes())
    render["sources"][0] = dict(path="/historical/worktree/" + MAIN,
                                **module.identity(new.read_bytes()))
    legacy.save_json(repo / module.RENDER, render)
    render_ref = dict(path="/historical/worktree/" + module.RENDER,
                      **module.identity((repo / module.RENDER).read_bytes()))
    for name in (module.PRINT, module.REVIEW):
        receipt = json.loads((repo / name).read_bytes())
        receipt["checked_records"] = [render_ref if row["path"].endswith("/" + module.RENDER)
                                      else row for row in receipt["checked_records"]]
        legacy.save_json(repo / name, receipt)
    # Retain the ledger revision actually bound by the render.
    (repo / module.LEDGER).write_text("Historical ledger\n")
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "explicit manuscript"], check=True)
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    return module, repo, revision, output


@pytest.mark.parametrize("explicit_stages", [False, True])
def test_new_main_exact_committed_bytes_and_relocated_cli(selected, tmp_path, explicit_stages):
    module, repo, revision, output = selected
    (repo / MAIN).write_text("Dirty manuscript must not enter export\n")
    stages = module.stage_paths() if explicit_stages else None
    result = module.build(repo, revision, revision, revision, output,
                          stages=stages, main_text=MAIN)
    index = json.loads((output / "REVIEW_INDEX.json").read_bytes())
    assert index["schema"] == "publication_direct_review_v3"
    assert index["main_text"] == MAIN
    assert index["stages"] == module.stage_paths()
    assert (output / MAIN).read_text().startswith("# Test")
    assert not (output / module.MAIN).exists()
    assert result["entrypoints"]["markdown"] == MAIN
    assert result["publication_ready"] is False
    relocated = tmp_path / "relocated-final-main"
    shutil.move(output, relocated)
    shutil.rmtree(repo)
    completed = subprocess.run([sys.executable, "-I", "-S", "-B",
        str(relocated / module.RUNNER), "verify", str(relocated),
        "--manifest-sha256", result["manifest"]["sha256"]], cwd=tmp_path,
        env={"PATH": "/no-git-or-original-checkout"}, capture_output=True,
        text=True, timeout=30)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == result


@pytest.mark.parametrize("name", [None, "", "../main.md", "/main.md", "a//main.md",
                                  "main.json", "benchmark_tools/results/PUBLICATION_PROGRESS.md",
                                  "benchmark_tools/PUBLICATION_REVIEW_COMPONENT.md", "LICENSE.md", 42])
def test_invalid_explicit_main(name):
    with pytest.raises(ValueError):
        legacy.bundle.manuscript_path(name)


def test_wrong_selected_main_rejected_without_output(selected):
    module, repo, revision, output = selected
    with pytest.raises(ValueError, match="Render source"):
        module.build(repo, revision, revision, revision, output,
                     main_text="benchmark_tools/results/different.md")
    assert not output.exists()


@pytest.mark.parametrize("defect", ["missing_main", "different_main", "legacy_v1", "legacy_v2",
                                    "missing_stages", "entrypoint", "missing_bound_render"])
def test_explicit_main_chain_cannot_be_reinterpreted(selected, defect):
    module, repo, revision, output = selected
    module.build(repo, revision, revision, revision, output, main_text=MAIN)
    path = output / "REVIEW_INDEX.json"
    index = json.loads(path.read_bytes())
    if defect == "missing_main":
        del index["main_text"]
    elif defect == "different_main":
        index["main_text"] = "benchmark_tools/results/different.md"
    elif defect.startswith("legacy_"):
        index["schema"] = "publication_direct_review_" + defect.split("_", 1)[1]
    elif defect == "missing_stages":
        del index["stages"]
    elif defect == "entrypoint":
        index["entrypoints"]["markdown"] = module.MAIN
    else:
        printed = output / module.PRINT
        data = json.loads(printed.read_bytes())
        data["checked_records"] = [row for row in data["checked_records"]
                                   if not row["path"].endswith("/" + module.RENDER)]
        legacy.save_json(printed, data)
        for row in index["files"]:
            if row["path"] == module.PRINT:
                row.update(module.identity(printed.read_bytes()))
    legacy.save_json(path, index)
    with pytest.raises(ValueError):
        module.verify(output, module.identity(path.read_bytes())["sha256"])


def test_build_cli_accepts_explicit_main(selected, tmp_path):
    module, repo, revision, output = selected
    completed = subprocess.run([sys.executable, "-I", "-S", "-B", module.__file__,
        "build", "--repo", str(repo), "--review-revision", revision,
        "--ledger-revision", revision, "--workflow-revision", revision,
        "--output", str(output), "--main-text", MAIN], cwd=tmp_path,
        capture_output=True, text=True, timeout=30)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["entrypoints"]["markdown"] == MAIN
