import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import bundle_publication_review as bundle


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


@pytest.fixture
def prepared(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    for key, value in (("user.name", "Test"), ("user.email", "test@example.invalid")):
        subprocess.run(["git", "-C", str(repo), "config", key, value], check=True)
    root = "/historical/worktree"
    html = "benchmark_tools/results/review.html"
    pdf = "benchmark_tools/results/printed/document.pdf"
    page = "benchmark_tools/results/review/page_001.png"
    asset = "benchmark_tools/results/figure.pdf"
    contents = {bundle.MAIN: b"# Test\n[Figure](figure.pdf)\n[Ledger](PUBLICATION_PROGRESS.md)\n",
                bundle.LEDGER: b"Historical ledger\n", asset: b"figure fixture, not a real PDF",
                html: b'<html><a href="figure.pdf">Figure</a><a href="PUBLICATION_PROGRESS.md">Ledger</a><a href="#refs">Refs</a><a href="https://example.invalid/">Citation</a></html>',
                pdf: b"PDF fixture, not a real PDF", page: b"page image fixture",
                "LICENSE.md": b"Test license", bundle.GUIDE: b"Test guide",
                bundle.RUNNER: Path(bundle.__file__).read_bytes()}
    for name, content in contents.items():
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    def record(name):
        return dict(path=root + "/" + name, **bundle.identity((repo / name).read_bytes()))

    render = dict(status="manuscript_review_rendered", publication_ready=False,
                  html=record(html), sources=[record(bundle.MAIN), dict(path="/usr/bin/pandoc", bytes=1, sha256="a" * 64)],
                  targets=[record(asset), record(bundle.LEDGER)], unique_targets=2,
                  occurrences=[dict(path=asset, url="figure.pdf"), dict(path=bundle.LEDGER, url="PUBLICATION_PROGRESS.md")],
                  local_occurrences=2)
    save_json(repo / bundle.RENDER, render)
    printed = dict(status="verified_html_printed", publication_ready=False, returncode=0,
                   page_count=1, pdf=record(pdf), checked_records=[record(bundle.RENDER), record(html)])
    save_json(repo / bundle.PRINT, printed)
    reviewed = dict(status="pdf_bounds_checked", publication_ready=False, page_count=1,
                    bounds_violations=[], checked_records=[record(pdf), record(bundle.RENDER)], rendered_pages=[record(page)])
    save_json(repo / bundle.REVIEW, reviewed)
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "retained review"], check=True)
    older = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    (repo / bundle.LEDGER).write_text("Updated ledger after rendering\n")
    subprocess.run(["git", "-C", str(repo), "add", bundle.LEDGER], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "later status"], check=True)
    newer = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    return repo, older, newer, tmp_path / "output"


def build(prepared):
    repo, older, newer, output = prepared
    result = bundle.build(repo, newer, older, newer, output)
    return output, result


def test_exports_exact_historical_ledger_and_ignores_dirty_files(prepared):
    repo = prepared[0]
    (repo / bundle.MAIN).write_text("Uncommitted unrelated draft\n")
    (repo / "untracked.txt").write_text("Must not enter component")
    output, result = build(prepared)
    assert (output / bundle.LEDGER).read_text() == "Historical ledger\n"
    assert (output / bundle.MAIN).read_text().startswith("# Test")
    assert result["status"] == "publication_direct_review_verified"
    assert result["direct_targets"] == 2
    assert result["local_html_occurrences"] == 2
    assert result["publication_ready"] is False
    assert result["transitive_evidence_included"] is False
    assert not (output / "usr/bin/pandoc").exists()


def test_wrong_ledger_revision_fails_without_creating_output(prepared):
    repo, _, newer, output = prepared
    with pytest.raises(ValueError, match="historical receipt"):
        bundle.build(repo, newer, newer, newer, output)
    assert not output.exists()


def test_offline_verifier_after_removing_source_repository(prepared, tmp_path):
    output, result = build(prepared)
    shutil.rmtree(prepared[0])
    relocated = tmp_path / "relocated"
    shutil.move(output, relocated)
    completed = subprocess.run([sys.executable, "-I", "-B", str(relocated / bundle.RUNNER),
                               "verify", str(relocated), "--manifest-sha256", result["manifest"]["sha256"]],
                              capture_output=True, text=True, check=True, env={"PATH": "/no-git"})
    assert json.loads(completed.stdout) == result


def test_repeat_output_refused(prepared):
    build(prepared)
    with pytest.raises(FileExistsError):
        build(prepared)


def custom_stages(prepared):
    repo, _, newer, output = prepared
    stages = {key: f"benchmark_tools/results/new_{key}.json" for key in ("render", "print", "review")}
    render = json.loads((repo / bundle.RENDER).read_text())
    save_json(repo / stages["render"], render)
    replacement = dict(path="/historical/worktree/" + stages["render"],
                       **bundle.identity((repo / stages["render"]).read_bytes()))
    for role, old in (("print", bundle.PRINT), ("review", bundle.REVIEW)):
        data = json.loads((repo / old).read_text())
        data["checked_records"] = [replacement if row["path"].endswith("/" + bundle.RENDER) else row
                                   for row in data["checked_records"]]
        save_json(repo / stages[role], data)
    subprocess.run(["git", "-C", str(repo), "add", *stages.values()], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "new selected stages"], check=True)
    latest = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    return repo, latest, prepared[1], newer, output, stages


def test_explicit_stages_and_offline_v2(prepared, tmp_path):
    repo, latest, ledger, workflow, output, stages = custom_stages(prepared)
    result = bundle.build(repo, latest, ledger, workflow, output, stages=stages)
    index = json.loads((output / "REVIEW_INDEX.json").read_text())
    assert index["schema"] == "publication_direct_review_v2"
    assert index["stages"] == stages
    assert not (output / bundle.RENDER).exists()
    assert (output / stages["render"]).exists()
    shutil.rmtree(repo)
    completed = subprocess.run([sys.executable, "-I", "-B", str(output / bundle.RUNNER),
        "verify", str(output), "--manifest-sha256", result["manifest"]["sha256"]],
        capture_output=True, text=True, check=True, env={"PATH": "/no-git"})
    assert json.loads(completed.stdout) == result


@pytest.mark.parametrize("stages", [
    {}, {"render": "a.json"},
    dict(render="../a.json", print="b.json", review="c.json"),
    dict(render="a.json", print="a.json", review="c.json"),
    dict(render="a.html", print="b.json", review="c.json"),
    dict(render="REVIEW_INDEX.json", print="b.json", review="c.json"),
    dict(render="a.json", print="b.json", review="c.json", extra="d.json")])
def test_bad_stage_paths_fail_before_export(prepared, stages):
    repo, older, newer, output = prepared
    with pytest.raises(ValueError):
        bundle.build(repo, newer, older, newer, output, stages=stages)
    assert not output.exists()


@pytest.mark.parametrize("kind", ["missing", "null", "different", "legacy_override"])
def test_stage_manifest_tampering(prepared, kind):
    repo, latest, ledger, workflow, output, stages = custom_stages(prepared)
    bundle.build(repo, latest, ledger, workflow, output, stages=stages)
    path = output / "REVIEW_INDEX.json"
    index = json.loads(path.read_text())
    if kind == "missing":
        del index["stages"]
    elif kind == "null":
        index["stages"] = None
    elif kind == "different":
        index["stages"]["render"] = "wrong.json"
    else:
        index["schema"] = "publication_direct_review_v1"
    save_json(path, index)
    with pytest.raises(ValueError):
        bundle.verify(output, bundle.identity(path.read_bytes())["sha256"])


def test_selected_render_chain_must_be_checked(prepared):
    repo, latest, ledger, workflow, output, stages = custom_stages(prepared)
    path = repo / stages["print"]
    data = json.loads(path.read_text())
    data["checked_records"] = [row for row in data["checked_records"]
                               if not row["path"].endswith("/" + stages["render"])]
    save_json(path, data)
    subprocess.run(["git", "-C", str(repo), "add", stages["print"]], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-qm", "unbound print"], check=True)
    latest = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    with pytest.raises(ValueError, match="Selected render/print/review chain"):
        bundle.build(repo, latest, ledger, workflow, output, stages=stages)


@pytest.mark.parametrize("kind", ["changed", "missing", "extra", "symlink", "parent_symlink", "mode"])
def test_payload_changes_rejected(prepared, tmp_path, kind):
    output, result = build(prepared)
    asset = output / "benchmark_tools/results/figure.pdf"
    if kind == "changed":
        asset.write_bytes(b"different")
    elif kind == "missing":
        asset.unlink()
    elif kind == "extra":
        (output / "extra.txt").write_text("extra")
    elif kind == "symlink":
        asset.unlink()
        asset.symlink_to(output / bundle.MAIN)
    elif kind == "parent_symlink":
        parent = asset.parent
        target = tmp_path / "outside"
        shutil.move(parent, target)
        parent.symlink_to(target, target_is_directory=True)
    else:
        asset.chmod(0o600)
    with pytest.raises(ValueError):
        bundle.verify(output, result["manifest"]["sha256"])


def test_coordinated_index_payload_tampering_needs_external_anchor(prepared):
    output, result = build(prepared)
    index = output / "REVIEW_INDEX.json"
    data = json.loads(index.read_text())
    asset = output / "benchmark_tools/results/figure.pdf"
    asset.write_bytes(b"changed")
    for row in data["files"]:
        if row["path"] == "benchmark_tools/results/figure.pdf":
            row.update(bundle.identity(asset.read_bytes()))
    save_json(index, data)
    with pytest.raises(ValueError, match="external digest"):
        bundle.verify(output, result["manifest"]["sha256"])
    with pytest.raises(ValueError, match="Historical review identity"):
        bundle.verify(output, bundle.identity(index.read_bytes())["sha256"])


@pytest.mark.parametrize("field", ["publication_ready", "redistribution_clearance", "transitive_evidence_included"])
def test_scope_cannot_be_upgraded(prepared, field):
    output, _ = build(prepared)
    index = output / "REVIEW_INDEX.json"
    data = json.loads(index.read_text())
    data[field] = True
    save_json(index, data)
    with pytest.raises(ValueError, match="scope"):
        bundle.verify(output, bundle.identity(index.read_bytes())["sha256"])


@pytest.mark.parametrize("url", ["/absolute.pdf", "../../../../outside", "%2Fabsolute.pdf", "file:///old/path", "javascript:alert(1)", "//remote.invalid/asset", "bad%5Cpath.pdf"])
def test_nonportable_or_escaping_links_rejected(url):
    with pytest.raises(ValueError):
        bundle.direct_links(('<a href="' + url + '">').encode(), bundle.MAIN)


@pytest.mark.parametrize("name", ["", "../outside", "/absolute", "a//b", "a/./b", "a\\b", None])
def test_invalid_relative_paths(name):
    with pytest.raises(ValueError):
        bundle.relative(name)


@pytest.mark.parametrize("change", ["failed_print", "bounds", "pages", "target_count", "occurrences", "conflict", "mapping"])
def test_inconsistent_stage_receipts_rejected(prepared, change):
    repo = prepared[0]
    render, printed, reviewed = [json.loads((repo / name).read_text()) for name in (bundle.RENDER, bundle.PRINT, bundle.REVIEW)]
    if change == "failed_print":
        printed["returncode"] = 1
    elif change == "bounds":
        reviewed["bounds_violations"] = [dict(page=1)]
    elif change == "pages":
        reviewed["page_count"] = 2
    elif change == "target_count":
        render["unique_targets"] = 3
    elif change == "occurrences":
        render["occurrences"] = []
    elif change == "mapping":
        render["occurrences"][0]["url"] = "wrong.pdf"
    else:
        reviewed["checked_records"][0]["bytes"] += 1
    with pytest.raises(ValueError):
        bundle.evidence(render, printed, reviewed)


def test_rendered_link_inventory_is_rechecked(prepared):
    output, result = build(prepared)
    index = output / "REVIEW_INDEX.json"
    data = json.loads(index.read_text())
    # A self-consistent replacement HTML still must agree with the target inventory.
    html_name = result["entrypoints"]["html"]
    content = (output / html_name).read_bytes().replace(b"figure.pdf", b"wrong.pdf")
    (output / html_name).write_bytes(content)
    replacement = bundle.identity(content)
    for row in data["files"]:
        if row["path"] == html_name:
            row.update(replacement)
    for name in (bundle.RENDER, bundle.PRINT):
        value = json.loads((output / name).read_text())
        if name == bundle.RENDER:
            value["html"].update(replacement)
        else:
            for row in value["checked_records"]:
                if row["path"].endswith("/" + html_name):
                    row.update(replacement)
                elif row["path"].endswith("/" + bundle.RENDER):
                    row.update(bundle.identity((output / bundle.RENDER).read_bytes()))
        save_json(output / name, value)
        for row in data["files"]:
            if row["path"] == name:
                row.update(bundle.identity((output / name).read_bytes()))
    value = json.loads((output / bundle.REVIEW).read_text())
    value["checked_records"][1].update(bundle.identity((output / bundle.RENDER).read_bytes()))
    save_json(output / bundle.REVIEW, value)
    for row in data["files"]:
        if row["path"] == bundle.REVIEW:
            row.update(bundle.identity((output / bundle.REVIEW).read_bytes()))
    save_json(index, data)
    with pytest.raises(ValueError, match="Rendered direct links"):
        bundle.verify(output, bundle.identity(index.read_bytes())["sha256"])
