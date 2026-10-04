"""Separate retained family evaluations from mentions and unproven tuning history."""

from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import inventory_development_families as inventory


def catalogs():
    names = ["RefOG%03d.txt" % i for i in range(1, 71)]
    ob = {"refogs": dict.fromkeys(names, {}), "development": names[::2], "validation": names[1::2]}
    swiss = {"families": ["S%02d" % i for i in range(18)]}
    return ob, swiss


def refog(name="RefOG001.txt"):
    return {"refog": name, "genes": 3, "true_positive": 1, "false_positive": 2, "false_negative": 2}


def swiss(name="S00"):
    return {"family": name, "counts_without_prior": {"TP": 1, "FP": 2, "FN": 3, "TN": 4}}


def test_nested_scores_aliases_and_mentions_have_distinct_scopes():
    aliases, roles = inventory.catalog(*catalogs())
    doc = {"partition": "validation", "all/variants~": {"score": {"refog_records": [refog("RefOG002")]}},
           "cells": [{"families": [swiss()]}], "families": ["S00", "RefOG002.txt"],
           "reference": {"S00": {"source": "not scored"}}, "other": "RefOG003.txt"}
    before = deepcopy(doc)
    blocks, mentions, unresolved = inventory.discover(doc, aliases)
    assert doc == before
    assert len(blocks) == 2
    assert blocks[0]["pointer"] == "/all~1variants~0/score/refog_records"
    assert blocks[0]["records"][0]["family"] == "RefOG002.txt"
    assert roles["RefOG002.txt"] == "validation"
    assert {b["dataset"] for b in blocks} == {"OrthoBench", "QfO_SwissTrees"}
    assert len(mentions) == 3
    assert unresolved == []
    assert not any(x["family"] == "RefOG003.txt" for x in mentions)


def test_no_scores_inferred_from_plans_statistics_or_family_count():
    aliases, _ = inventory.catalog(*catalogs())
    blocks, mentions, unresolved = inventory.discover({"planned": {"families": ["S00", "S01"]},
        "families": 70, "reference_families": 70, "annotation": {"S00": {}},
        "contrast": [{"family": "S00", "F1": .2}]}, aliases)
    assert blocks == [] and len(mentions) == 3
    assert unresolved == []


def test_subset_ordinal_candidate_label_is_not_canonical_family_evidence():
    aliases, _ = inventory.catalog(*catalogs())
    diagnostic = {"audit": {"refog_records": [{"refog": "RefOG001", "genes": 19,
                  "possible_pairs": 171, "recovered_pairs": 21}]}}
    blocks, mentions, unresolved = inventory.discover(diagnostic, aliases)
    assert not blocks and not mentions
    assert unresolved == [{"pointer": "/audit/refog_records", "kind": "ordinal_candidate_diagnostic",
                           "recorded_labels": ["RefOG001"], "reference_family_mapping_established": False}]
    with pytest.raises(ValueError): inventory.discover({"refog_records": [{"refog": "RefOG001"}]}, aliases)


def test_synthetic_fixture_labels_are_retained_without_biological_namespace():
    aliases, _ = inventory.catalog(*catalogs())
    blocks, mentions, unresolved = inventory.discover({"refog_records": [refog("family0.txt")]}, aliases)
    assert not blocks and not mentions
    assert unresolved[0]["kind"] == "unmapped_scored_reference_labels"
    assert unresolved[0]["recorded_labels"] == ["family0.txt"]
    assert unresolved[0]["reference_family_mapping_established"] is False


@pytest.mark.parametrize("damage", ["missing", "overlap", "wrong_size", "duplicate_swiss", "empty_swiss"])
def test_catalog_invariants(damage):
    ob, sw = catalogs()
    if damage == "missing": ob["refogs"].pop("RefOG070.txt")
    elif damage == "overlap": ob["validation"][0] = ob["development"][0]
    elif damage == "wrong_size": ob["development"].append(ob["validation"].pop())
    elif damage == "duplicate_swiss": sw["families"][1] = sw["families"][0]
    elif damage == "empty_swiss": sw["families"][0] = ""
    with pytest.raises(ValueError): inventory.catalog(ob, sw)


@pytest.mark.parametrize("namespace", ["OrthoBench", "QfO_SwissTrees"])
@pytest.mark.parametrize("damage", ["unknown", "duplicate", "negative", "nan", "boolean", "empty"])
def test_malformed_explicit_counts_cannot_silently_enter_inventory(namespace, damage):
    aliases, _ = inventory.catalog(*catalogs())
    row = refog() if namespace == "OrthoBench" else swiss()
    records = [row]
    if damage == "unknown": row["refog" if namespace == "OrthoBench" else "family"] = "wrong"
    elif damage == "duplicate": records.append(deepcopy(row))
    elif damage == "empty": records = []
    else:
        target, key = (row, "true_positive") if namespace == "OrthoBench" else (row["counts_without_prior"], "TP")
        target[key] = -1 if damage == "negative" else float("nan") if damage == "nan" else True
    with pytest.raises(ValueError): inventory.scored_records(records, namespace, aliases)


def git(repo, *argv):
    return subprocess.check_output(["git", *argv], cwd=repo, text=True).strip()


def test_actual_git_snapshot_ignores_dirty_files_and_retains_duplicate_evidence(tmp_path, monkeypatch):
    git(tmp_path, "init", "-q")
    (tmp_path / "benchmark_tools/results").mkdir(parents=True)
    ob, sw = catalogs()
    paths = {"OrthoBench": "benchmark_tools/results/ob.json", "QfO_SwissTrees": "benchmark_tools/results/sw.json"}
    for name, doc in (("OrthoBench", ob), ("QfO_SwissTrees", sw)):
        (tmp_path / paths[name]).write_text(json.dumps(doc))
    evidence = {"partition": "validation", "score": {"refog_records": [refog("RefOG002.txt")]},
                "cells": [{"families": [swiss()]}]}
    for name in ("first", "duplicate"):
        (tmp_path / ("benchmark_tools/results/" + name + ".json")).write_text(json.dumps(evidence))
    git(tmp_path, "add", "benchmark_tools")
    git(tmp_path, "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid",
        "commit", "-qm", "Frozen results")
    revision = git(tmp_path, "rev-parse", "HEAD")
    monkeypatch.setattr(inventory, "CATALOGS", {name: (path, hashlib.sha256((tmp_path / path).read_bytes()).hexdigest())
                                               for name, path in paths.items()})
    # No mutable result read may replace the committed evidence.
    (tmp_path / paths["OrthoBench"]).write_text("not JSON")
    (tmp_path / "benchmark_tools/results/untracked.json").write_text("not JSON")
    report = inventory.write(tmp_path, tmp_path / "output", revision)
    assert report["snapshot_commit"] == revision and report["json_files_scanned"] == 4
    assert report["summary"]["OrthoBench"]["scored_blocks"] == 2
    assert report["summary"]["OrthoBench"]["distinct_count_vectors"] == 1
    assert report["summary"]["QfO_SwissTrees"]["scored_blocks"] == 2
    assert len(inventory.association_rows(report)) == 4
    validation = next(r for r in report["family_rows"] if r["family"] == "RefOG002.txt")
    assert validation["original_partition"] == "validation" and validation["publication_exposure"] == "development_exposed"
    assert validation["independent_native_experiments"] is None
    assert validation["causal_tuning_influence_established"] is False
    assert report["family_disjoint_validation_established"] is False
    for name in ("inventory.json", "inventory.md", "family_evidence.tsv"):
        assert (tmp_path / "output" / name).is_file()
    with pytest.raises(FileExistsError): inventory.write(tmp_path, tmp_path / "output", revision)
    # Freeze additional retained local metadata explicitly; do not silently read it.
    (tmp_path / "benchmark_tools/results/untracked.json").write_text(json.dumps(evidence))
    manifest = tmp_path / "local_manifest.json"
    local = inventory.freeze_local_reports(tmp_path, manifest, revision)
    assert len(local["local_reports"]) == 1
    manifest_sha = hashlib.sha256(manifest.read_bytes()).hexdigest()
    expanded = inventory.collect(tmp_path, revision, manifest, manifest_sha)
    assert expanded["json_files_scanned"] == 5 and expanded["local_json_bytes"] > 0
    assert expanded["summary"]["OrthoBench"]["scored_blocks"] == 3
    with pytest.raises(ValueError): inventory.collect(tmp_path, revision, manifest, "wrong")
    (tmp_path / "benchmark_tools/results/untracked.json").write_text("changed")
    with pytest.raises(ValueError): inventory.collect(tmp_path, revision, manifest, manifest_sha)


def test_existing_and_broken_symlink_output_refused(tmp_path):
    with pytest.raises(FileExistsError): inventory.write(tmp_path, tmp_path)
    target = tmp_path / "broken"
    target.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError): inventory.write(tmp_path, target)


def test_actual_selected_inventory():
    repo = Path(__file__).resolve().parents[2]
    output = repo / "benchmark_tools/results/development_family_inventory_20261004"
    if not (output / "inventory.json").exists():
        pytest.skip("Prepare source before selected snapshot collection")
    report = json.loads((output / "inventory.json").read_text())
    assert report["snapshot_commit"] == inventory.SNAPSHOT
    assert report["json_files_scanned"] == 1899
    assert report["snapshot_json_bytes"] == 144384742
    assert report["local_json_bytes"] == 108265645
    assert len(report["family_rows"]) == 88
    assert report["summary"]["OrthoBench"]["families_with_scored_evidence"] == 70
    assert report["summary"]["QfO_SwissTrees"]["families_with_scored_evidence"] == 18
    assert all(r["publication_exposure"] == "development_exposed" for r in report["family_rows"])
    assert len({r["family"] for r in report["family_rows"] if r["original_partition"] == "validation"}) == 35
    with (output / "family_evidence.tsv").open(newline="") as stream:
        actual = list(csv.DictReader(stream, delimiter="\t"))
    rows = inventory.association_rows(report)
    assert actual == [{k: "NA" if v is None else str(v) for k, v in row.items()} for row in rows]
