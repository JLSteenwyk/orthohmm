import copy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import audit_qfo_factorial_swiss as retained_audit
from benchmark_tools import bind_native_qfo_swiss_uncertainty as module
from benchmark_tools import bootstrap_qfo_factorial as bootstrap
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_audit_qfo_factorial_swiss import synthetic
from tests.unit.test_audit_native_qfo_swiss_counts import fixture as ordinary_fixture
from tests.unit.test_audit_recovered_native_qfo_swiss_counts import fixture as recovered_fixture, run
from tests.unit.test_export_native_qfo_factorial_scores import write


@pytest.fixture(scope="module")
def retained_bootstrap(tmp_path_factory):
    root = tmp_path_factory.mktemp("synthetic_retained_bootstrap")
    entries, baseline = synthetic(root)
    counts = retained_audit.assemble(entries, baseline)
    result = bootstrap.bootstrap(counts)
    result["status"] = "paired_corrected_qfo_factorial_swiss_intervals"
    return counts, result


def fixture(tmp_path, monkeypatch, retained_bootstrap, different_counts=False):
    refs = recovered_fixture(tmp_path, monkeypatch, different_counts=different_counts)
    audit = run(refs)
    audit_ref = write(tmp_path / "audit.json", audit)
    counts = json.loads(Path(refs[1]["path"]).read_text())
    frozen = copy.deepcopy(retained_bootstrap[1])
    provenance = tmp_path / "synthetic_provenance"
    provenance.write_text("synthetic test provenance\n")
    frozen.update(counts=refs[1], source=record(provenance), protocol=record(provenance),
        corrected_protocol=record(provenance), helpers=[])
    frozen_ref = write(tmp_path / "bootstrap.json", frozen)
    monkeypatch.setattr(module, "BOOTSTRAP_SHA", frozen_ref["sha256"])
    return refs[0], audit_ref, refs[1], frozen_ref, counts, audit


def bind(refs):
    snapshot, audit, counts, frozen, *_ = refs
    return module.bind(snapshot["path"], snapshot["sha256"], [(audit["path"], audit["sha256"])],
                       counts["path"], frozen["path"])


def native_rows(counts, indexes):
    rows = {}
    for index in indexes:
        row = copy.deepcopy(counts["cells"][index])
        row.update(index=index + 6, native_job_id=10, count_audit={}, admission={},
            native_endpoint_f1=row["aggregate"]["F1"], retained_family_records_identical=True,
            retained_aggregate_identical=True)
        rows[row["cell"]] = row
    return rows


def test_only_matched_two_cell_contrast_reuses_original_intervals(retained_bootstrap):
    counts, frozen = retained_bootstrap
    rows = native_rows(counts, [0, 1])
    bindings, contrasts = module.project(rows, counts, frozen)
    assert len(bindings) == 2 and len(contrasts) == 14
    matched = [row for row in contrasts if row["status"] == "native_records_matched"]
    assert len(matched) == 1 and matched[0]["name"] == "R_at_P0_C0"
    original = next(row for row in frozen["comparisons"] if row["name"] == "R_at_P0_C0")
    assert matched[0]["metrics"] == original["metrics"]
    assert matched[0]["family_differences"] == original["family_differences"]
    assert all(row["metrics"] is None for row in contrasts if row["status"] != "native_records_matched")


def test_aggregate_agreement_does_not_establish_family_match(retained_bootstrap):
    counts, frozen = retained_bootstrap
    rows = native_rows(counts, [0, 1])
    row = rows["p0_c0_r1"]
    row["families"][0]["represented_genes"][0] = "different_member"
    row["retained_family_records_identical"] = False
    bindings, contrasts = module.project(rows, counts, frozen)
    assert bindings["p0_c0_r1"]["status"] == "native_records_differ"
    assert not any(c["status"] == "native_records_matched" for c in contrasts)
    effect = next(c for c in contrasts if c["name"] == "R_at_P0_C0")
    assert effect["differing_cells"] == ["p0_c0_r1"] and effect["metrics"] is None


def test_false_equality_claim_rejected(retained_bootstrap):
    counts, frozen = retained_bootstrap
    rows = native_rows(counts, [0])
    rows["p0_c0_r0"]["aggregate"]["F1"] = .123
    with pytest.raises(ValueError, match="equality flags"):
        module.project(rows, counts, frozen)


@pytest.mark.parametrize("field,value", [("multiplicity_endpoints", 3), ("seed", 1),
    ("alpha", .1), ("replicates", 100), ("quantile_method", "nearest"), ("publication_ready", True)])
def test_frozen_bootstrap_scope_cannot_change(retained_bootstrap, field, value):
    counts, frozen = retained_bootstrap
    changed = copy.deepcopy(frozen)
    changed[field] = value
    with pytest.raises(ValueError):
        module.project(native_rows(counts, [0, 1]), counts, changed)


@pytest.mark.parametrize("change", ["weights", "order", "metrics", "points", "difference"])
def test_retained_definitions_and_arithmetic_required(retained_bootstrap, change):
    counts, frozen = retained_bootstrap
    changed = copy.deepcopy(frozen)
    if change == "weights":
        changed["comparisons"][0]["weights"][0] = 0
    elif change == "order":
        changed["comparisons"].reverse()
    elif change == "metrics":
        changed["comparisons"][0]["metrics"].pop("TPR")
    elif change == "points":
        changed["point_estimates"]["p0_c0_r0"]["F1"] = .123
    else:
        changed["comparisons"][0]["metrics"]["F1"]["difference"] = .123
    with pytest.raises(ValueError):
        module.project(native_rows(counts, [0, 1]), counts, changed)


def test_recovered_single_cell_binding_adds_no_contrast_or_timing(tmp_path, monkeypatch, retained_bootstrap):
    result = bind(fixture(tmp_path, monkeypatch, retained_bootstrap))
    assert list(result["bound_cells"]) == ["p0_c0_r1"]
    assert result["replicates_reused"] == 100000 and result["multiplicity_endpoints"] == 42
    assert result["new_bootstrap_draws"] == 0
    assert not any(row["status"] == "native_records_matched" for row in result["contrasts"])
    for key in ("publication_ready", "new_accuracy_or_resource_admission", "independent_confirmation"):
        assert result[key] is False


def test_actual_differing_audited_counts_are_not_forced_to_cached_values(tmp_path, monkeypatch, retained_bootstrap):
    result = bind(fixture(tmp_path, monkeypatch, retained_bootstrap, different_counts=True))
    assert result["bound_cells"]["p0_c0_r1"]["status"] == "native_records_differ"
    assert all(row["metrics"] is None for row in result["contrasts"])


def test_ordinary_count_audit_can_bind_without_recount(tmp_path, monkeypatch, retained_bootstrap):
    snapshot_ref, retained_ref, _ = ordinary_fixture(tmp_path, monkeypatch)
    count_audit = module.ordinary.audit(snapshot_ref["path"], snapshot_ref["sha256"], retained_ref["path"])
    audit_ref = write(tmp_path / "ordinary_audit.json", count_audit)
    snapshot = json.loads(Path(snapshot_ref["path"]).read_text())
    admission = snapshot["rows"][0]["admission"]
    combined = tmp_path / "scientific_snapshot"
    module.reporter.export(snapshot["plan"]["path"], snapshot["plan"]["sha256"],
        [(admission["path"], admission["sha256"])], [], combined)
    combined_ref = record(combined / "report.json")
    frozen = copy.deepcopy(retained_bootstrap[1])
    provenance = tmp_path / "synthetic_provenance"
    provenance.write_text("synthetic test provenance\n")
    frozen.update(counts=retained_ref, source=record(provenance), protocol=record(provenance),
        corrected_protocol=record(provenance), helpers=[])
    frozen_ref = write(tmp_path / "bootstrap.json", frozen)
    monkeypatch.setattr(module, "BOOTSTRAP_SHA", frozen_ref["sha256"])
    result = module.bind(combined_ref["path"], combined_ref["sha256"],
        [(audit_ref["path"], audit_ref["sha256"])], retained_ref["path"], frozen_ref["path"])
    assert result["bound_cells"]["p0_c0_r0"]["status"] == "native_records_matched"
    assert all(row["metrics"] is None for row in result["contrasts"])


@pytest.mark.parametrize("field,value", [("schema", "cached_count_audit"), ("source", {}),
    ("publication_ready", True), ("historical_intervals_attached", True), ("new_bootstrap_draws", 1),
    ("independent_confirmation", True), ("new_accuracy_or_resource_admission", True),
    ("reference_relation_count", 0), ("count_conversion", "raw counts without prior")])
def test_changed_count_scope_rejected(tmp_path, monkeypatch, retained_bootstrap, field, value):
    refs = list(fixture(tmp_path, monkeypatch, retained_bootstrap))
    refs[5][field] = value
    refs[1] = write(Path(refs[1]["path"]), refs[5])
    with pytest.raises(ValueError):
        bind(refs)


@pytest.mark.parametrize("field,value", [("admission", {}), ("index", 6), ("native_job_id", 11),
    ("native_endpoint_f1", .123), ("resources", {}), ("timing_admitted", True), ("timing_eligible", True)])
def test_admission_binding_and_failed_timing_preserved(tmp_path, monkeypatch, retained_bootstrap, field, value):
    refs = list(fixture(tmp_path, monkeypatch, retained_bootstrap))
    refs[5]["cells"][0][field] = value
    refs[1] = write(Path(refs[1]["path"]), refs[5])
    with pytest.raises(ValueError):
        bind(refs)


def test_duplicate_counts_rejected(tmp_path, monkeypatch, retained_bootstrap):
    refs = fixture(tmp_path, monkeypatch, retained_bootstrap)
    with pytest.raises(ValueError, match="Duplicate"):
        module.bind(refs[0]["path"], refs[0]["sha256"], [(refs[1]["path"], refs[1]["sha256"])] * 2,
                    refs[2]["path"], refs[3]["path"])


def test_native_raw_binding_is_checked_without_recount(tmp_path, monkeypatch, retained_bootstrap):
    refs = fixture(tmp_path, monkeypatch, retained_bootstrap)
    Path(refs[5]["cells"][0]["raw_file"]["path"]).write_bytes(b"changed")
    with pytest.raises(ValueError):
        bind(refs)


def test_no_overwrite_before_binding(tmp_path, monkeypatch):
    output = tmp_path / "existing"
    output.write_text("retain me")
    monkeypatch.setattr(sys, "argv", ["bind", "--snapshot", "absent", "--snapshot-sha256", "unused",
        "--counts-audit", "absent", "unused", "--retained-counts", "absent", "--bootstrap", "absent",
        "--output", str(output)])
    with pytest.raises(ValueError, match="Output already exists"):
        module.main()
    assert output.read_text() == "retain me"
