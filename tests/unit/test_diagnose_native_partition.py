"""Invented coverage fixtures; no selected benchmark data or inference."""

import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import diagnose_native_partition as diagnostic


def save(path, value):
    path.write_text(json.dumps(value))
    bindings = diagnostic.Bindings()
    bindings.bind(path)
    return bindings.files[str(path)]


def fixture(tmp_path):
    root = tmp_path / "run"
    inputs = root / "input"
    originals = tmp_path / "originals"
    working = root / "native/orthohmm_working_res"
    phylogeny = root / "native/orthohmm_phylogeny"
    checkpoint = working / "high_sensitivity_checkpoint"
    for path in (inputs, originals, checkpoint, phylogeny):
        path.mkdir(parents=True)
    pins = []
    for filename, contents in (("alpha.fa", ">a description\nAAAA\n>c\nCCCC\n"),
                               ("beta.fa", ">b\nBBBB\n")):
        (inputs / filename).write_text(contents)
        (originals / filename).write_text(contents)
        bindings = diagnostic.Bindings()
        bindings.bind(originals / filename)
        pins.append(bindings.files[str(originals / filename)])
    order = ["beta.fa", "alpha.fa"]
    digest = hashlib.sha256()
    for filename, gene in (("beta.fa", "b"), ("alpha.fa", "a"), ("alpha.fa", "c")):
        digest.update(json.dumps([filename, gene], separators=(",", ":")).encode() + b"\n")
    save(root / "preparation.json", dict(genes=3, gene_ownership_sha256=digest.hexdigest(),
                                        per_species_counts={"alpha.fa": 2, "beta.fa": 1}))
    names = checkpoint / "gene_names.txt"
    names.write_text("a\nb\nc\n")
    bindings = diagnostic.Bindings()
    bindings.bind(names)
    save(checkpoint / "manifest.json", dict(genes=3, files={"gene_names.txt": {
        k: bindings.files[str(names)][k] for k in ("bytes", "sha256")}}))
    (phylogeny / "orthohmm_root_hogs.tsv").write_text(
        "root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\ta,c\n"
        "RootHOG0000001\tFamily0000001\tb\n")
    (working / "orthohmm_edges_clustered.txt").write_text("a c\nb\n")
    (root / "native/orthohmm_orthogroups.txt").write_text("OG0: a c\nOG1: b\n")
    save(phylogeny / "provenance_manifest.json", dict(input_cluster_sha256=
         hashlib.sha256(b"a c\nb\n").hexdigest()))
    save(phylogeny / "reconciliation_summary.json", dict(candidate_families=2, root_hogs=2))
    run = dict(index=0, output_root=str(root), inputs=pins, native_order=order, genes=3, proteomes=2)
    plan = save(tmp_path / "plan.json", dict(runs=[run]))
    failure = save(tmp_path / "failure.json", dict(index=0, job_id=999, plan=plan,
                   status="terminal_factorial_review_failed", terminal_reviewed=False,
                   next_identity_authorized=False, accuracy_evaluated=False, automatic_retry=False))
    return root, plan, failure


def test_full_invented_diagnostic_preserves_failure_without_admission(tmp_path):
    _, plan, failure = fixture(tmp_path)
    result = diagnostic.diagnose(plan, 0, failure)
    assert result["input_genes"] == 3
    assert result["input_proteomes"] == 2
    assert result["frozen_root_coverage_gate"] == dict(status="passed", groups=2)
    assert result["parsers_identical"] is True
    assert result["all_output_partitions_identical"] is True
    assert result["source_payload_matches_native_input"] is True
    assert result["source_family_ids_complete"] is True
    assert result["root_count_matches_summary"] is True
    assert result["checkpoint_names_lexical"] is True
    assert all(stage["complete_unique_input_partition"] for stage in result["stages"].values())
    assert all(result[key] is False for key in (
        "accuracy_evaluated", "terminal_reviewed", "next_identity_authorized",
        "native_outputs_validated", "automatic_retry", "resources_admitted",
        "historical_failure_cause_established"))
    bindings = diagnostic.Bindings()
    for ref in result["evidence"]:
        bindings.bind(ref["path"], ref)


@pytest.mark.parametrize("genes,missing,extra,duplicates", [
    (["a", "b"], ["c"], [], {}),
    (["a", "b", "c", "x"], [], ["x"], {}),
    (["a", "a", "b", "c"], [], [], {"a": 2}),
    ([], ["a", "b", "c"], [], {}),
])
def test_missing_extra_duplicate_are_distinct(genes, missing, extra, duplicates):
    report = diagnostic.coverage({"g": genes}, {"a": "s1", "b": "s2", "c": "s1"})
    assert report["missing"] == missing
    assert report["extra"] == extra
    assert report["duplicates"] == duplicates
    assert report["memberships"] == len(genes)
    assert report["complete_unique_input_partition"] is False


def test_per_species_memberships_do_not_hide_duplicates():
    report = diagnostic.coverage({"x": ["a", "a", "b", "extra"]}, {"a": "s1", "b": "s2", "c": "s1"})
    assert report["per_species"]["s1"] == dict(input=2, observed_unique=1, memberships=2, missing=1)
    assert report["per_species"]["s2"] == dict(input=1, observed_unique=1, memberships=1, missing=0)


@pytest.mark.parametrize("text", [
    "bad\theader\n", "root_hog\tsource_family\tgenes\ng\tf\n",
    "root_hog\tsource_family\tgenes\ng\tf\ta\ng\tf\tb\n",
    "root_hog\tsource_family\tgenes\ng\t\ta\n",
    "root_hog\tsource_family\tgenes\ng\tf\ta,\n",
    "root_hog\tsource_family\tgenes\ng\tf\t a\n",
])
def test_structurally_invalid_roots_fail(tmp_path, text):
    path = tmp_path / "roots.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        diagnostic.root_groups(path)


@pytest.mark.parametrize("text", ["no colon\n", ": a\n", "OG: \n", "OG: a\nOG: b\n"])
def test_invalid_named_groups_fail(tmp_path, text):
    path = tmp_path / "groups.txt"
    path.write_text(text)
    with pytest.raises(ValueError):
        diagnostic.text_groups(path, named=True)


def test_partition_digest_ignores_names_and_order_not_multiplicity():
    first = diagnostic.canonical_partition_digest({"a": ["c", "b"], "d": ["a"]})
    assert first == diagnostic.canonical_partition_digest({"x": ["a"], "y": ["b", "c"]})
    assert first != diagnostic.canonical_partition_digest({"x": ["a"], "y": ["b", "b", "c"]})
    assert first != diagnostic.canonical_partition_digest({"x": ["a", "b", "c"]})


@pytest.mark.parametrize("key", ["plan", "index", "terminal_reviewed", "next_identity_authorized",
                                  "automatic_retry", "accuracy_evaluated", "status"])
def test_failure_identity_and_disposition_cannot_be_changed(tmp_path, key):
    _, plan, failure = fixture(tmp_path)
    value = json.loads(Path(failure["path"]).read_text())
    value[key] = {"plan": {}, "index": 1, "status": "success"}.get(key, True)
    failure = save(Path(failure["path"]), value)
    with pytest.raises(ValueError):
        diagnostic.diagnose(plan, 0, failure)


def test_missing_root_gene_is_reported_and_frozen_gate_rejects(tmp_path):
    root, plan, failure = fixture(tmp_path)
    path = root / "native/orthohmm_phylogeny/orthohmm_root_hogs.tsv"
    path.write_text(path.read_text().replace("a,c", "a"))
    report = diagnostic.diagnose(plan, 0, failure)
    assert report["stages"]["root_hogs"]["missing"] == ["c"]
    assert report["frozen_root_coverage_gate"]["status"] == "failed"
    assert report["all_output_partitions_identical"] is False
    assert report["source_payload_matches_native_input"] is False
    assert report["native_outputs_validated"] is False


@pytest.mark.parametrize("mutation", ["input_bytes", "extra_file", "checkpoint_hash", "preparation_counts"])
def test_provenance_mismatch_cannot_be_diagnosed_as_valid_input(tmp_path, mutation):
    root, plan, failure = fixture(tmp_path)
    if mutation == "input_bytes":
        (root / "input/alpha.fa").write_text(">a\nTTTT\n>c\nCCCC\n")
    elif mutation == "extra_file":
        (root / "input/unexpected.fa").write_text(">x\nXXXX\n")
    elif mutation == "checkpoint_hash":
        (root / "native/orthohmm_working_res/high_sensitivity_checkpoint/gene_names.txt").write_text("a\nb\n")
    else:
        path = root / "preparation.json"
        value = json.loads(path.read_text())
        value["per_species_counts"]["alpha.fa"] += 1
        save(path, value)
    with pytest.raises(ValueError):
        diagnostic.diagnose(plan, 0, failure)


def test_binding_rejects_later_mutation(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("original")
    bindings = diagnostic.Bindings()
    bindings.bind(path)
    path.write_text("mutated")
    with pytest.raises(ValueError, match="Pinned evidence differs"):
        bindings.finish()


def test_binding_rejects_symlink(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("original")
    link = tmp_path / "linked.txt"
    link.symlink_to(path)
    with pytest.raises(ValueError, match="direct absolute"):
        diagnostic.Bindings().bind(link)
