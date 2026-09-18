import gzip
from io import StringIO
import json
import subprocess
import sys

import pytest

from benchmark_tools import prepare_qfo_corrected_orthofinder_pairs as module
from benchmark_tools.audit_orthofinder_pair_tables import audit_tables
from benchmark_tools.orthofinder_to_pairwise import iter_pairs
from tests.unit.test_audit_orthofinder_pair_tables import fixture


def test_native_relations_not_mcl_cliques(tmp_path):
    owners, member = fixture(tmp_path)
    groups = {"OG0": list(owners)}
    content = {"results_directory": str(tmp_path), "native_pairs": audit_tables(tmp_path, owners, "ABC", member)}
    native, clique = StringIO(), StringIO()
    assert module.write_native(content, owners, groups, native) == (2, 0)
    assert module.write_groups(groups, owners, clique) == (5, 0)
    native_pairs = {tuple(row.split("\t")) for row in native.getvalue().splitlines()}
    assert native_pairs == set(iter_pairs(tmp_path))
    assert "b\tc" not in native.getvalue().splitlines()
    assert "b\tc" in clique.getvalue().splitlines()


def test_mcl_equivalent_to_existing_group_converter(tmp_path):
    owners, _ = fixture(tmp_path)
    groups = {"OG0": list(owners)}
    for species in "ABC":
        (tmp_path / f"{species}.fasta").write_text("".join(
            f">{gene}\nMKA\n" for gene, owner in owners.items() if owner == species))
    partition = tmp_path / "groups.txt"
    partition.write_text(" ".join(owners) + "\n")
    native_script = module.Path(module.__file__).resolve().parents[1] / "qfo_benchmark/og_to_pairwise.py"
    old = subprocess.check_output([sys.executable, str(native_script), str(partition), str(tmp_path)], text=True)
    new = StringIO()
    module.write_groups(groups, owners, new)
    assert set(old.splitlines()) == set(new.getvalue().splitlines())
    assert len(new.getvalue().splitlines()) == len(set(new.getvalue().splitlines()))


def test_native_duplicate_removal_matches_existing_converter_set(tmp_path):
    owners, member = fixture(tmp_path)
    path = tmp_path / "Orthologues/Orthologues_A/A__v__B.tsv"
    with path.open("a") as stream:
        stream.write("OG0\tsp|a|A\tsp|b|B\n")
    content = {"results_directory": str(tmp_path), "native_pairs": audit_tables(tmp_path, owners, "ABC", member)}
    stream = StringIO()
    assert module.write_native(content, owners, {"OG0": list(owners)}, stream) == (2, 1)
    assert {tuple(s.split("\t")) for s in stream.getvalue().splitlines()} == set(iter_pairs(tmp_path))


@pytest.mark.parametrize("kind", ["counts", "order", "total", "duplicates"])
def test_native_binding_changes_rejected(tmp_path, kind):
    owners, member = fixture(tmp_path)
    content = {"results_directory": str(tmp_path), "native_pairs": audit_tables(tmp_path, owners, "ABC", member)}
    audit = content["native_pairs"]
    if kind == "counts":
        audit["species_pairs"][0]["forward"]["rows"] += 1
    elif kind == "order":
        audit["species_pairs"].reverse()
    elif kind == "total":
        audit["distinct_pairs"] += 1
    else:
        audit["converter_duplicate_pairs"] += 1
    with pytest.raises(ValueError):
        module.write_native(content, owners, {"OG0": list(owners)}, StringIO())


@pytest.mark.parametrize("kind", ["missing", "duplicate", "alias", "unknown"])
def test_invalid_partition_rejected(tmp_path, kind):
    owners, _ = fixture(tmp_path)
    groups = {"OG0": list(owners)}
    if kind == "missing":
        groups["OG0"].pop()
    elif kind == "duplicate":
        groups["extra"] = [groups["OG0"][0]]
    elif kind == "alias":
        owners["tr|a|other"] = "A"
        groups["OG0"].append("tr|a|other")
    else:
        groups["OG0"].append("unknown")
    with pytest.raises(ValueError):
        module.write_groups(groups, owners, StringIO())


def admission():
    return {"status": "corrected_orthofinder_native_evidence_admitted", "accuracy_evaluated": False,
        "publication_ready": False, "scheduler": {"State": "COMPLETED", "ExitCode": "0:0",
        "AllocCPUS": "32", "NodeList": "bizon"}, "content": {"genes": 984137, "species": 78,
        "native_pairs": {"status": "native_orthofinder_pair_tables_verified", "mcl_membership_checked": True,
        "directed_tables": 6006, "species": 78, "input_genes": 984137, "distinct_pairs": 10,
        "converter_emitted_pairs": 12, "converter_duplicate_pairs": 2}}}


def test_corrected_admission_scope():
    value = admission()
    assert module.validate_admission(value) == value["content"]


@pytest.mark.parametrize("key,value", [("genes", 976504), ("species", 77)])
def test_original_or_incomplete_input_refused(key, value):
    report = admission()
    report["content"][key] = value
    with pytest.raises(ValueError):
        module.validate_admission(report)


@pytest.mark.parametrize("key,value", [("mcl_membership_checked", False), ("directed_tables", 3003),
    ("distinct_pairs", True), ("distinct_pairs", 0), ("converter_duplicate_pairs", 1)])
def test_incomplete_pair_audit_refused(key, value):
    report = admission()
    report["content"]["native_pairs"][key] = value
    with pytest.raises(ValueError):
        module.validate_admission(report)


def test_pending_admission_refused_before_io(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_JOB_ID", "200")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n101|PENDING|0:0|0:00|bizon|2\n")
    with pytest.raises(ValueError, match="COMPLETED"):
        module.prepare(tmp_path, module.METHODS[0], tmp_path / "missing", "0" * 64, 101)


def prepared_fixture(tmp_path, monkeypatch, valid_ids=("a", "a2", "b", "c")):
    # Mock only upstream admission/runtime identity; execute writers, hashing and filtering.
    owners, member = fixture(tmp_path)
    native = audit_tables(tmp_path, owners, "ABC", member)
    groups = {"OG0": list(owners)}
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    inputs = tmp_path / "input"
    inputs.mkdir()
    for species in "ABC":
        (inputs / f"{species}.fasta").write_text("".join(
            f">{g}\nMKA\n" for g, s in owners.items() if s == species))
    executor = tmp_path / "benchmarks/work/publication_qfo_corrected_orthofinder_admission_v1/benchmark_tools"
    executor.mkdir(parents=True)
    source = executor / "admit_qfo_corrected_orthofinder.py"
    source.write_text("# fixture admission\n")
    plan = results / "qfo_corrected_primary_commands_20260918.json"
    plan.write_text(json.dumps({"inputs": [module.record(p) for p in sorted(inputs.iterdir())],
                               "input_directory": str(inputs)}))
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": {g: i for i, g in enumerate(valid_ids)}}, stream)
    env = results / "qfo_assessment_environment_20260917.json"
    env.write_text(json.dumps({"reference_files": [module.record(mapping)]}))
    content = {"native_pairs": native, "results_directory": str(tmp_path), "checkpoint_groups": 1,
               "checkpoint": {"path": "fixture-checkpoint"}, "sequence_ids": {"path": "fixture-ids"}}
    report = {"source": module.record(source), "content": content,
              "checked_records": [module.record(plan), *native["checked_files"]]}
    path = tmp_path / "admission.json"
    path.write_text(json.dumps(report))
    monkeypatch.setenv("SLURM_JOB_ID", "200")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setattr(module, "PLAN_SHA", module.record(plan)["sha256"])
    monkeypatch.setattr(module, "ENV_SHA", module.record(env)["sha256"])
    monkeypatch.setattr(module, "validate_admission", lambda report: report["content"])
    monkeypatch.setattr(module, "input_universe", lambda dataset: (owners, list("ABC")))
    monkeypatch.setattr(module, "read_checkpoint", lambda *args: groups)
    monkeypatch.setattr(module.subprocess, "check_output", lambda argv, **kw:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n101|COMPLETED|0:0|0:01|bizon|2\n"
        if argv[0] == "sacct" else module.ADMITTER + "\n")
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: None)
    return path, module.record(path)["sha256"]


@pytest.mark.parametrize("method,expected", [(module.METHODS[0], 2), (module.METHODS[1], 5)])
def test_conversion_lifecycle_and_no_overwrite(tmp_path, monkeypatch, method, expected):
    path, sha = prepared_fixture(tmp_path, monkeypatch)
    result = module.prepare(tmp_path, method, path, sha, 101)
    assert result["status"] == "corrected_orthofinder_pairs_prepared_unscored"
    assert result["total_pairs"] == result["retained_pairs"] == expected
    assert result["semantics"] == module.SEMANTICS[method]
    assert result["pairs"]["sha256"] == result["filtered_pairs"]["sha256"]
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, method, path, sha, 101)


def test_mapping_loss_preserves_failure_evidence(tmp_path, monkeypatch):
    path, sha = prepared_fixture(tmp_path, monkeypatch, valid_ids=("a", "a2"))
    with pytest.raises(ValueError, match="mapping loss"):
        module.prepare(tmp_path, module.METHODS[0], path, sha, 101)
    output = tmp_path / "benchmarks/results/qfo_corrected_comparator_pairs_v1/orthofinder_full"
    assert json.loads((output / "results.json").read_text())["status"] == "failed"
    assert (output / "pairs.partial.tsv").exists()
    assert not (output / "pairs.tsv").exists()


def test_native_mutation_during_conversion_refused(tmp_path, monkeypatch):
    path, sha = prepared_fixture(tmp_path, monkeypatch)
    original = module.filter_pairs

    def mutate(*args):
        value = original(*args)
        (tmp_path / "Orthologues/Orthologues_A/A__v__B.tsv").write_text("changed\n")
        return value

    monkeypatch.setattr(module, "filter_pairs", mutate)
    with pytest.raises(ValueError):
        module.prepare(tmp_path, module.METHODS[0], path, sha, 101)
    output = tmp_path / "benchmarks/results/qfo_corrected_comparator_pairs_v1/orthofinder_full"
    assert json.loads((output / "results.json").read_text())["status"] == "failed"
    assert not (output / "pairs.tsv").exists()
