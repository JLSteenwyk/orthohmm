import json
import subprocess
import sys

import pytest

from benchmark_tools import audit_orthomcl_failure_membership as module
from benchmark_tools.audit_orthomcl_native_groups import audit as group_audit
from benchmark_tools.audit_orthomcl_reference_impact import failed_query_targets
from tests.unit.test_audit_orthomcl_native_groups import fixture


def search():
    return dict(input_proteins=7, failed_queries=3, failed_queries_with_query_hits=0,
        failed_queries_with_subject_hits=2, diagnostics=[dict(gene=g, length=10,
            query_failed=g != "E", has_query_hits=g == "E", has_subject_hits=g in ("A", "G", "E"),
            has_self_hit=g == "E", messages=[dict(category="setup_failure" if g != "E" else "selenocysteine_replacement")])
            for g in ("A", "D", "G", "E")])


def test_grouped_singleton_unindexed_and_incoming_are_distinct(tmp_path):
    result = module.describe(search(), *fixture(tmp_path))
    rows = {r["gene"]: r for r in result["records"]}
    assert result["failed_queries_by_membership"] == dict(final_group=1, mcl_singleton=1, absent_from_index=1)
    assert rows["A"]["final_group_line"] == 1 and rows["A"]["final_group_size"] == 3
    assert rows["D"]["native_cluster_id"] == 1 and rows["D"]["final_group_line"] is None
    assert rows["G"]["native_cluster_id"] is None and rows["G"]["has_subject_hits"] is True
    assert rows["E"]["query_failed"] is False and rows["E"]["final_group_size"] == 2
    failed, targets = failed_query_targets(result, {"A": 1, "D": 2, "G": 3}, allow_grouped=True)
    assert len(failed) == 3 and targets == {1, 2, 3}
    with pytest.raises(ValueError, match="grouped exposure"):
        failed_query_targets(result, {"A": 1, "D": 2, "G": 3})


@pytest.mark.parametrize("problem", ["universe", "unknown", "duplicate", "bool", "length", "outgoing", "self", "category", "total", "incoming"])
def test_bad_search_evidence_rejected(tmp_path, problem):
    value = search()
    row = value["diagnostics"][0]
    if problem == "universe":
        value["input_proteins"] += 1
    elif problem == "unknown":
        row["gene"] = "unknown"
    elif problem == "duplicate":
        value["diagnostics"].append(row)
    elif problem == "bool":
        row["query_failed"] = 1
    elif problem == "length":
        row["length"] = 0
    elif problem == "outgoing":
        row["has_query_hits"] = True
    elif problem == "self":
        row["has_self_hit"] = True
    elif problem == "category":
        row["messages"] = []
    else:
        value["failed_queries" if problem == "total" else "failed_queries_with_subject_hits"] += 1
    with pytest.raises(ValueError):
        module.describe(value, *fixture(tmp_path))


def admitted_fixture(tmp_path):
    paths = [p.rename(tmp_path / n) for p, n in zip(fixture(tmp_path),
             ("all_orthomcl.out", "all_ortho.mcl", "all_ortho.idx", "all.gg"))]
    group_path = tmp_path / "groups.json"
    result = group_audit(*paths, group_path)
    native = dict(status="recovered_orthomcl_native_outputs_admitted", conversion_authorized=True,
        query_coverage=search(), outputs=[module.record(group_path)], content=result["content"],
        native_groups=module.record(paths[0]))
    report = dict(status="recovered_orthomcl_search_evidence_verified", search_admitted=True, query_coverage=search())
    return native, report


@pytest.mark.parametrize("problem", [None, "status", "coverage", "group_pin", "content", "ambiguous"])
def test_pinned_admission_flow(tmp_path, problem):
    native, report = admitted_fixture(tmp_path)
    if problem == "status":
        native["conversion_authorized"] = False
    elif problem == "coverage":
        native["query_coverage"]["failed_queries"] += 1
    elif problem == "group_pin":
        native["native_groups"]["sha256"] = "0" * 64
    elif problem == "content":
        native["content"]["final_groups"] += 1
    elif problem == "ambiguous":
        native["outputs"] *= 2
    paths = [tmp_path / n for n in ("search.json", "native.json")]
    for path, value in zip(paths, (report, native)):
        path.write_text(json.dumps(value))
    args = (paths[0], module.record(paths[0])["sha256"], paths[1], module.record(paths[1])["sha256"])
    if problem:
        with pytest.raises(ValueError):
            module.audit(*args)
    else:
        result = module.audit(*args)
        assert result["failed_queries"] == 3
        assert result["accuracy_admitted"] is False and result["publication_ready"] is False
        assert result["checked_records"]


def test_cli_writes_content_audit_once(tmp_path):
    native, report = admitted_fixture(tmp_path)
    paths = [tmp_path / n for n in ("search.json", "native.json")]
    for path, value in zip(paths, (report, native)):
        path.write_text(json.dumps(value))
    output = tmp_path / "exposure.json"
    command = [sys.executable, module.__file__, "--search", str(paths[0]), "--native", str(paths[1]),
               "--search-sha256", module.record(paths[0])["sha256"],
               "--native-sha256", module.record(paths[1])["sha256"], "--output", str(output)]
    subprocess.run(command, check=True, capture_output=True)
    result = json.loads(output.read_text())
    assert result["status"] == "recovered_failure_membership_content_verified"
    before = module.record(output)
    assert subprocess.run(command, capture_output=True).returncode != 0
    assert module.record(output) == before


@pytest.mark.parametrize("problem", [None, "no_opt_in", "verification", "changed", "record"])
def test_native_representation_requires_verified_opt_in(tmp_path, monkeypatch, problem):
    native, report = admitted_fixture(tmp_path)
    report.update(status="recovered_search_native_representation_verified",
                  database_representation=dict(exact_sequence_parity=False,
                      limitations=["Native O deletions retained; impact not measured."]))
    paths = [tmp_path / n for n in ("search.json", "native.json")]
    for path, value in zip(paths, (report, native)):
        path.write_text(json.dumps(value))
    evidence = module.record(paths[0])
    if problem == "record":
        evidence["sha256"] = "0" * 64
    calls = []

    def verify(root, path, digest):
        calls.append((root, path, digest))
        if problem == "verification":
            raise ValueError("Representation verification failed")
        verified = dict(report)
        if problem == "changed":
            verified["search_admitted"] = False
        return verified, [], [evidence], {}, ""

    monkeypatch.setattr(module, "verify_admission", verify)
    args = (paths[0], module.record(paths[0])["sha256"], paths[1], module.record(paths[1])["sha256"])
    root = None if problem == "no_opt_in" else tmp_path
    if problem:
        with pytest.raises(ValueError):
            module.audit(*args, native_representation_root=root)
    else:
        result = module.audit(*args, native_representation_root=root)
        assert result["database_representation"]["exact_sequence_parity"] is False
        assert report["database_representation"]["limitations"][0] in result["limitations"]
        assert result["failed_queries"] == 3
        assert result["accuracy_admitted"] is False
        assert evidence in result["checked_records"]
    assert calls == ([] if root is None else [(root, paths[0], args[1])])
