"""Synthetic admission handoffs with real Swiss counting and interval kernels."""

from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import audit_composed_native_qfo_swiss_counts as audit
from benchmark_tools import bind_composed_native_qfo_swiss_uncertainty as binder
from benchmark_tools import bootstrap_qfo_factorial as bootstrap
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_audit_native_qfo_swiss_counts import fixture as ordinary_fixture
from tests.unit.test_bind_native_qfo_swiss_uncertainty import retained_bootstrap
from tests.unit.test_export_native_qfo_factorial_scores import write


def fixture(tmp_path, monkeypatch, problem=None, matching=True):
    snapshot_ref, counts_ref, raw_ref = ordinary_fixture(tmp_path, monkeypatch,
        problem=problem, different_counts=not matching)
    counts = json.loads(Path(counts_ref["path"]).read_text())
    if matching:
        for key in ("families", "aggregate"):
            counts["cells"][7][key] = deepcopy(counts["cells"][0][key])
        counts["cells"][7]["raw_file"] = deepcopy(counts["cells"][0]["raw_file"])
        counts_ref = write(Path(counts_ref["path"]), counts)
        monkeypatch.setattr(audit.ordinary, "RETAINED_COUNTS_SHA", counts_ref["sha256"])
    snapshot = json.loads(Path(snapshot_ref["path"]).read_text())
    row = deepcopy(snapshot["rows"][0])
    row.update(index=12, cell="p1_c1_r1", status="supplied_composed_native_admission", accuracy_admitted=True)
    snapshot["rows"][0].update(status="no_supplied_native_admission", admission=None, accuracy_admitted=False)
    snapshot["rows"][6] = row
    for item in snapshot["rows"]:
        item.setdefault("accuracy_admitted", False)
    snapshot.update(schema="composed_native_qfo_scientific_reporting_v1",
        source=record(audit.reporter.__file__), supplied_admissions=1, final_admission=row["admission"])
    replay = {k: v for k, v in snapshot.items() if k not in ("source", "outputs")}
    # Exercise the new interfaces/count kernel without pretending these invented
    # admissions pass the separately tested production resource/runtime route.
    monkeypatch.setattr(audit.reporter, "collect", lambda *args: deepcopy(replay))
    snapshot_ref = write(Path(snapshot_ref["path"]), snapshot)
    return snapshot_ref, counts_ref, raw_ref, snapshot


def run(refs, index=12):
    return audit.audit(refs[0]["path"], refs[0]["sha256"], refs[1]["path"], index)


def test_selected_raw_uses_original_count_kernel_and_preserves_reference(tmp_path, monkeypatch):
    refs = fixture(tmp_path, monkeypatch)
    result = run(refs)
    assert result["schema"] == "composed_native_qfo_swiss_family_count_audit_v1"
    assert result["selected_index"] == 12 and len(result["cells"]) == 1
    assert len(result["families"]) == 18 and result["reference_relation_count"] == 270
    cell = result["cells"][0]
    assert cell["cell"] == "p1_c1_r1" and cell["raw_file"] == refs[2]
    assert cell["retained_family_records_identical"] and cell["retained_aggregate_identical"]
    assert cell["differing_families"] == []
    assert result["new_bootstrap_draws"] == 0 and not result["historical_intervals_attached"]
    assert not result["new_accuracy_or_resource_admission"] and not result["publication_ready"]


def test_real_count_difference_not_replaced_with_retained_values(tmp_path, monkeypatch):
    cell = run(fixture(tmp_path, monkeypatch, matching=False))["cells"][0]
    assert not cell["retained_family_records_identical"]
    assert len(cell["differing_families"]) == 18


@pytest.mark.parametrize("problem", ["truth", "members", "duplicate", "inventory", "family_metric",
    "missing_metric", "missing_raw", "ambiguous_raw", "unchecked_execution", "unchecked_raw"])
def test_unchecked_or_inconsistent_raw_refused(tmp_path, monkeypatch, problem):
    with pytest.raises(ValueError):
        run(fixture(tmp_path, monkeypatch, problem=problem))


@pytest.mark.parametrize("index", [True, 9, 11, 13, "12"])
def test_only_explicit_integer_composed_selection(index):
    with pytest.raises(ValueError, match="composed native12 index"):
        audit.audit("absent", "unused", "absent", index)


@pytest.mark.parametrize("field,value", [("schema", "native_qfo_reporting_snapshot_v1"),
    ("source", {}), ("supplied_admissions", 7), ("publication_ready", True)])
def test_resealed_snapshot_does_not_bypass_report_replay(tmp_path, monkeypatch, field, value):
    refs = fixture(tmp_path, monkeypatch)
    changed = deepcopy(refs[3])
    changed[field] = value
    with pytest.raises(ValueError):
        run((write(Path(refs[0]["path"]), changed), *refs[1:]))


@pytest.mark.parametrize("target", ["raw", "table", "retained"])
def test_post_export_bound_file_tampering_rejected(tmp_path, monkeypatch, target):
    refs = fixture(tmp_path, monkeypatch)
    path = Path(refs[{"raw": 2, "retained": 1}.get(target, 0)]["path"])
    if target == "table":
        path = path.with_name("scores.tsv")
    path.write_bytes(b"changed")
    with pytest.raises(ValueError):
        run(refs)


def test_selection_does_not_recount_other_admitted_cells(tmp_path, monkeypatch):
    refs = fixture(tmp_path, monkeypatch)
    calls = []
    original = audit.ordinary.family_cell

    def observed(*args):
        calls.append(args[0])
        return original(*args)

    monkeypatch.setattr(audit.ordinary, "family_cell", observed)
    run(refs)
    assert calls == [refs[2]]


def test_unavailable_selected_cell_is_not_filled(tmp_path, monkeypatch):
    refs = fixture(tmp_path, monkeypatch)
    snapshot = deepcopy(refs[3])
    snapshot["rows"][6].update(status="no_supplied_native_admission", accuracy_admitted=False)
    replay = {k: v for k, v in snapshot.items() if k not in ("source", "outputs")}
    monkeypatch.setattr(audit.reporter, "collect", lambda *args: deepcopy(replay))
    refs = (write(Path(refs[0]["path"]), snapshot), *refs[1:])
    with pytest.raises(ValueError, match="no composed admission"):
        run(refs)


@pytest.fixture(scope="module")
def frozen(retained_bootstrap):
    counts = deepcopy(retained_bootstrap[0])
    for key in ("families", "aggregate"):
        counts["cells"][7][key] = deepcopy(counts["cells"][0][key])
    result = bootstrap.bootstrap(counts)
    result["status"] = "paired_corrected_qfo_factorial_swiss_intervals"
    return result


def binding_fixture(tmp_path, monkeypatch, frozen, matching=True):
    refs = fixture(tmp_path, monkeypatch, matching=matching)
    count_audit = run(refs)
    audit_ref = write(tmp_path / "audit.json", count_audit)
    provenance = tmp_path / "provenance"
    provenance.write_text("invented fixture provenance\n")
    value = deepcopy(frozen)
    value.update(counts=refs[1], source=record(provenance), protocol=record(provenance),
        corrected_protocol=record(provenance), helpers=[])
    frozen_ref = write(tmp_path / "bootstrap.json", value)
    monkeypatch.setattr(binder.original, "BOOTSTRAP_SHA", frozen_ref["sha256"])
    return refs, audit_ref, frozen_ref, count_audit


def bind(refs, duplicate=False):
    inputs, audit_ref, frozen_ref, _ = refs
    audits = [(audit_ref["path"], audit_ref["sha256"])] * (2 if duplicate else 1)
    return binder.bind(inputs[0]["path"], inputs[0]["sha256"], audits,
        inputs[1]["path"], frozen_ref["path"])


def test_composed_binding_keeps_full_adjustment_and_unavailable_contrasts(tmp_path, monkeypatch, frozen):
    result = bind(binding_fixture(tmp_path, monkeypatch, frozen))
    assert list(result["bound_cells"]) == ["p1_c1_r1"]
    assert result["bound_cells"]["p1_c1_r1"]["status"] == "native_records_matched"
    assert len(result["contrasts"]) == 14 and all(r["metrics"] is None for r in result["contrasts"])
    assert result["replicates_reused"] == 100000 and result["multiplicity_endpoints"] == 42
    assert not result["publication_ready"] and result["new_bootstrap_draws"] == 0


@pytest.mark.parametrize("field,value", [("schema", "cached"), ("source", {}), ("selected_index", True),
    ("selected_index", 11), ("new_bootstrap_draws", True), ("historical_intervals_attached", True),
    ("publication_ready", True), ("count_conversion", "changed"), ("reference_relation_count", 1)])
def test_resealed_count_scope_cannot_authorize_intervals(tmp_path, monkeypatch, frozen, field, value):
    refs = list(binding_fixture(tmp_path, monkeypatch, frozen))
    refs[3][field] = value
    refs[1] = write(Path(refs[1]["path"]), refs[3])
    with pytest.raises(ValueError):
        bind(refs)


@pytest.mark.parametrize("field,value", [("admission", {}), ("index", 11), ("native_job_id", 11),
    ("native_endpoint_f1", .123), ("raw_file", {}), ("retained_family_records_identical", False)])
def test_cell_binding_and_equality_claims_checked(tmp_path, monkeypatch, frozen, field, value):
    refs = list(binding_fixture(tmp_path, monkeypatch, frozen))
    refs[3]["cells"][0][field] = value
    refs[1] = write(Path(refs[1]["path"]), refs[3])
    with pytest.raises(ValueError):
        bind(refs)


def test_duplicate_audit_rejected(tmp_path, monkeypatch, frozen):
    with pytest.raises(ValueError, match="Duplicate"):
        bind(binding_fixture(tmp_path, monkeypatch, frozen), duplicate=True)


def test_binding_checks_raw_binding_without_recount(tmp_path, monkeypatch, frozen):
    refs = binding_fixture(tmp_path, monkeypatch, frozen)
    monkeypatch.setattr(audit.ordinary, "family_cell", lambda *args: pytest.fail("unexpected raw recount"))
    bind(refs)
    Path(refs[0][2]["path"]).write_bytes(b"changed")
    with pytest.raises(ValueError):
        bind(refs)


def mixed_fixture(tmp_path, monkeypatch, frozen, changed=False, matching=True):
    refs = binding_fixture(tmp_path, monkeypatch, frozen, matching=matching)
    inputs, composed_ref, frozen_ref, count_audit = refs
    counts = json.loads(Path(inputs[1]["path"]).read_text())
    old = deepcopy(count_audit)
    old.update(schema="allocated_native_qfo_swiss_family_count_audit_v1",
        status="supplied_allocated_native_swiss_family_counts_verified",
        selected_index=10, source=record(binder.allocated.allocated.__file__))
    row = deepcopy(counts["cells"][5])
    row.update(index=10, native_job_id=20, admission=inputs[3]["rows"][6]["admission"],
        native_endpoint_f1=row["aggregate"]["F1"], retained_family_records_identical=True,
        retained_aggregate_identical=True)
    old["cells"] = [row]
    old["checked_inputs"].append(row["raw_file"])
    old_ref = write(tmp_path / "allocated_audit.json", old)
    snapshot = deepcopy(inputs[3])
    snapshot["rows"][4].update(status="supplied_allocated_native_admission", accuracy_admitted=True,
        native_job_id=20, admission=row["admission"])
    snapshot["rows"][4]["scores"]["SwissTrees"] = row["native_endpoint_f1"]
    replay = {k: v for k, v in snapshot.items() if k not in ("source", "outputs")}
    monkeypatch.setattr(audit.reporter, "collect", lambda *args: deepcopy(replay))
    snapshot_ref = write(tmp_path / "mixed_snapshot.json", snapshot)
    if changed:
        count_audit["cells"][0]["families"][0]["represented_genes"][0] = "different_member"
        count_audit["cells"][0]["retained_family_records_identical"] = False
        composed_ref = write(Path(composed_ref["path"]), count_audit)
    return snapshot_ref, inputs[1], composed_ref, old_ref, frozen_ref, old


def mixed_bind(refs):
    snapshot, counts, composed_ref, old_ref, frozen_ref, _ = refs
    return binder.bind(snapshot["path"], snapshot["sha256"],
        [(r["path"], r["sha256"]) for r in (composed_ref, old_ref)], counts["path"], frozen_ref["path"])


def test_mixed_routes_expose_only_supported_candidate_effect(tmp_path, monkeypatch, frozen):
    result = mixed_bind(mixed_fixture(tmp_path, monkeypatch, frozen))
    matched = [r for r in result["contrasts"] if r["status"] == "native_records_matched"]
    assert len(matched) == 1 and matched[0]["name"] == "C_at_P1_R1"
    expected = next(r for r in frozen["comparisons"] if r["name"] == "C_at_P1_R1")
    assert matched[0]["metrics"] == expected["metrics"]
    assert matched[0]["family_differences"] == expected["family_differences"]
    assert result["multiplicity_endpoints"] == 42 and result["new_bootstrap_draws"] == 0
    assert "p1_c1_r0" not in result["bound_cells"]


def test_equal_aggregate_with_different_family_members_cannot_reuse_interval(tmp_path, monkeypatch, frozen):
    result = mixed_bind(mixed_fixture(tmp_path, monkeypatch, frozen, changed=True))
    assert result["bound_cells"]["p1_c1_r1"]["status"] == "native_records_differ"
    assert all(r["metrics"] is None for r in result["contrasts"])


def recovered_fixture(tmp_path, monkeypatch, frozen):
    refs = list(mixed_fixture(tmp_path, monkeypatch, frozen))
    snapshot = json.loads(Path(refs[0]["path"]).read_text())
    counts = json.loads(Path(refs[1]["path"]).read_text())
    old = deepcopy(refs[-1])
    old.update(schema="recovered_native_qfo_swiss_family_count_audit_v1",
        status="supplied_recovered_native_swiss_family_counts_verified",
        source=record(binder.original.recovered.__file__))
    row = deepcopy(counts["cells"][1])
    row.update(index=7, native_job_id=21, admission=snapshot["rows"][6]["admission"],
        native_endpoint_f1=row["aggregate"]["F1"], retained_family_records_identical=True,
        retained_aggregate_identical=True, resources=None, timing_admitted=False, timing_eligible=False)
    old["cells"] = [row]
    old["checked_inputs"].append(row["raw_file"])
    refs[3] = write(tmp_path / "recovered_audit.json", old)
    snapshot["rows"][1].update(status="supplied_recovered_scientific_admission", accuracy_admitted=True,
        native_job_id=21, admission=row["admission"])
    snapshot["rows"][1]["scores"]["SwissTrees"] = row["native_endpoint_f1"]
    replay = {k: v for k, v in snapshot.items() if k not in ("source", "outputs")}
    monkeypatch.setattr(audit.reporter, "collect", lambda *args: deepcopy(replay))
    refs[0] = write(Path(refs[0]["path"]), snapshot)
    refs[-1] = old
    return refs


def test_recovered_science_does_not_repair_failed_timing(tmp_path, monkeypatch, frozen):
    result = mixed_bind(recovered_fixture(tmp_path, monkeypatch, frozen))
    assert set(result["bound_cells"]) == {"p1_c1_r1", "p0_c0_r1"}
    assert all(r["metrics"] is None for r in result["contrasts"])


@pytest.mark.parametrize("field,value", [("replicates", 100), ("seed", 1), ("alpha", .1),
    ("multiplicity_endpoints", 3), ("quantile_method", "nearest"), ("publication_ready", True)])
def test_new_binding_cannot_reduce_or_change_retained_bootstrap_scope(tmp_path, monkeypatch, frozen, field, value):
    refs = list(binding_fixture(tmp_path, monkeypatch, frozen))
    value_ref = refs[2]
    changed = json.loads(Path(value_ref["path"]).read_text())
    changed[field] = value
    refs[2] = write(Path(value_ref["path"]), changed)
    monkeypatch.setattr(binder.original, "BOOTSTRAP_SHA", refs[2]["sha256"])
    with pytest.raises(ValueError, match="bootstrap scope"):
        bind(refs)


def test_ordinary_audit_route_remains_available_without_recount(tmp_path, monkeypatch, frozen):
    refs = list(mixed_fixture(tmp_path, monkeypatch, frozen))
    snapshot = json.loads(Path(refs[0]["path"]).read_text())
    counts = json.loads(Path(refs[1]["path"]).read_text())
    old = deepcopy(refs[-1])
    old.update(schema="native_qfo_swiss_family_count_audit_v1",
        status="supplied_native_swiss_family_counts_verified", source=record(binder.original.ordinary.__file__))
    row = deepcopy(counts["cells"][0])
    row.update(index=6, native_job_id=22, admission=snapshot["rows"][6]["admission"],
        native_endpoint_f1=row["aggregate"]["F1"], retained_family_records_identical=True,
        retained_aggregate_identical=True)
    old["cells"] = [row]
    old["checked_inputs"].append(row["raw_file"])
    refs[3] = write(tmp_path / "ordinary_audit.json", old)
    snapshot["rows"][0].update(status="supplied_native_admission", accuracy_admitted=True,
        native_job_id=22, admission=row["admission"])
    snapshot["rows"][0]["scores"]["SwissTrees"] = row["native_endpoint_f1"]
    replay = {k: v for k, v in snapshot.items() if k not in ("source", "outputs")}
    monkeypatch.setattr(audit.reporter, "collect", lambda *args: deepcopy(replay))
    monkeypatch.setattr(audit.ordinary, "family_cell", lambda *args: pytest.fail("unexpected raw recount"))
    refs[0] = write(Path(refs[0]["path"]), snapshot)
    result = mixed_bind(refs)
    assert set(result["bound_cells"]) == {"p1_c1_r1", "p0_c0_r0"}
    assert all(r["metrics"] is None for r in result["contrasts"])


@pytest.mark.parametrize("index,status", [(9, "no_supplied_native_admission"),
    (11, "retained_composed_scoring_failure")])
def test_missing_or_failed_other_cells_never_gain_counts(tmp_path, monkeypatch, frozen, index, status):
    refs = list(binding_fixture(tmp_path, monkeypatch, frozen))
    inputs = refs[0]
    snapshot = deepcopy(inputs[3])
    snapshot["rows"][index - 6].update(status=status, accuracy_admitted=False)
    replay = {k: v for k, v in snapshot.items() if k not in ("source", "outputs")}
    monkeypatch.setattr(audit.reporter, "collect", lambda *args: deepcopy(replay))
    snapshot_ref = write(Path(inputs[0]["path"]), snapshot)
    refs[0] = (snapshot_ref, *inputs[1:])
    cell = refs[3]["cells"][0]
    cell.update(index=index, cell=snapshot["rows"][index - 6]["cell"])
    refs[1] = write(Path(refs[1]["path"]), refs[3])
    with pytest.raises(ValueError):
        bind(refs)


@pytest.mark.parametrize("field,value", [("resources", {}), ("timing_admitted", True), ("timing_eligible", True)])
def test_mixed_route_never_repairs_recovered_failed_timing(tmp_path, monkeypatch, frozen, field, value):
    refs = recovered_fixture(tmp_path, monkeypatch, frozen)
    refs[-1]["cells"][0][field] = value
    refs[3] = write(Path(refs[3]["path"]), refs[-1])
    with pytest.raises(ValueError, match="failed timing"):
        mixed_bind(refs)


@pytest.mark.parametrize("module", [audit, binder])
def test_cli_refuses_existing_output_before_work(tmp_path, monkeypatch, module):
    output = tmp_path / "existing"
    output.write_text("retain me")
    extra = ["--index", "12"] if module is audit else ["--bootstrap", "absent", "--counts-audit", "absent", "unused"]
    monkeypatch.setattr(sys, "argv", ["test", "--snapshot", "absent", "--snapshot-sha256", "unused",
        "--retained-counts", "absent", "--output", str(output), *extra])
    with pytest.raises(ValueError, match="Output already exists"):
        module.main()
    assert output.read_text() == "retain me"
