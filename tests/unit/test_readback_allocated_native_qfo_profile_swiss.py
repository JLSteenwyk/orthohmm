"""Invented native-pair/rational readback fixtures, not biological truth."""

from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import readback_allocated_native_qfo_profile_swiss as reader
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_allocated_native_qfo_swiss_uncertainty import frozen, mixed_fixture, mixed_bind
from tests.unit.test_bind_native_qfo_swiss_uncertainty import retained_bootstrap
from tests.unit.test_export_native_qfo_factorial_scores import write


def fixture(tmp_path, monkeypatch, frozen):
    refs = mixed_fixture(tmp_path, monkeypatch, frozen)
    binding = mixed_bind(refs)
    binding_ref = write(tmp_path / "binding.json", binding)
    prior = deepcopy(binding)
    prior.update(schema="native_qfo_retained_swiss_uncertainty_binding_v1",
        source=record(Path(reader.__file__).with_name("bind_native_qfo_swiss_uncertainty.py")))
    for row in prior["contrasts"]:
        row.update(status="native_records_unavailable", metrics=None, family_differences=None)
    prior_ref = write(tmp_path / "prior.json", prior)
    monkeypatch.setattr(reader, "RELATIONS", 270)
    return refs[2], binding_ref, prior_ref


def run(refs):
    return reader.verify(*[v for r in refs for v in (r["path"], r["sha256"])])


def mutate(refs, target, change):
    refs = list(refs)
    value = json.loads(Path(refs[target]["path"]).read_text())
    change(value)
    refs[target] = write(Path(refs[target]["path"]), value)
    return refs


def test_independent_rational_raw_readback_and_interval_source(tmp_path, monkeypatch, frozen):
    result = run(fixture(tmp_path, monkeypatch, frozen))
    assert result["native_family_records_checked"] == 36 and result["families_checked"] == 18
    assert result["raw_rows_checked"] == 810 and result["profile_pair_labels_matched"] == 270
    assert result["contrast"]["name"] == "P_at_C0_R1"
    assert result["prior_matched_contrasts_checked"] == 0
    assert result["new_bootstrap_draws"] == 0 and not result["publication_ready"]
    assert not result["new_accuracy_or_resource_admission"]


@pytest.mark.parametrize("target,key,value", [(0, "schema", "cached"), (0, "source", {}),
    (0, "selected_index", True), (0, "historical_intervals_attached", True),
    (1, "schema", "native_qfo_retained_swiss_uncertainty_binding_v1"), (1, "source", {}),
    (1, "multiplicity_endpoints", 3), (1, "replicates_reused", 100), (1, "seed_reused", 1),
    (1, "alpha", .1), (1, "new_bootstrap_draws", True), (1, "publication_ready", True),
    (2, "schema", "cached"), (2, "source", {}), (2, "bootstrap", {})])
def test_resealed_scope_or_source_change_rejected(tmp_path, monkeypatch, frozen, target, key, value):
    refs = fixture(tmp_path, monkeypatch, frozen)
    with pytest.raises(ValueError):
        run(mutate(refs, target, lambda d: d.update({key: value})))


@pytest.mark.parametrize("change", ["difference", "interval", "weights", "family", "wins", "missing"])
def test_profile_contrast_cannot_be_resealed_to_new_values(tmp_path, monkeypatch, frozen, change):
    refs = fixture(tmp_path, monkeypatch, frozen)

    def alter(value):
        row = next(r for r in value["contrasts"] if r["name"] == "P_at_C0_R1")
        if change == "weights":
            row["weights"][0] = 1
        elif change == "family":
            row["family_differences"][0]["F1"] += .1
        elif change == "missing":
            row["missing_cells"] = ["p0_c0_r1"]
        elif change == "interval":
            row["metrics"]["F1"]["bonferroni_percentile_ci"][0] += .1
        elif change == "wins":
            row["metrics"]["F1"]["family_wins"] += 1
        else:
            row["metrics"]["F1"]["difference"] += .1

    with pytest.raises(ValueError):
        run(mutate(refs, 1, alter))


@pytest.mark.parametrize("target", [0, 1, 2])
def test_changed_supplied_digest_rejected(tmp_path, monkeypatch, frozen, target):
    refs = list(fixture(tmp_path, monkeypatch, frozen))
    refs[target] = dict(refs[target], sha256="0" * 64)
    with pytest.raises(ValueError, match="supplied report"):
        run(refs)


def test_raw_file_mutation_rejected_before_arithmetic(tmp_path, monkeypatch, frozen):
    refs = fixture(tmp_path, monkeypatch, frozen)
    audit = json.loads(Path(refs[0]["path"]).read_text())
    Path(audit["cells"][0]["raw_file"]["path"]).write_bytes(b"changed")
    with pytest.raises(ValueError):
        run(refs)


def test_prior_matched_interval_cannot_change(tmp_path, monkeypatch, frozen):
    refs = fixture(tmp_path, monkeypatch, frozen)

    def alter(value):
        row = value["contrasts"][0]
        row.update(status="native_records_matched", metrics={"F1": {"difference": 999}})

    with pytest.raises(ValueError, match="Prior matched contrast changed"):
        run(mutate(refs, 2, alter))


def test_no_overwrite_before_readback(tmp_path, monkeypatch):
    output = tmp_path / "existing"
    output.write_text("retain me")
    monkeypatch.setattr(sys, "argv", ["test", "--audit", "absent", "--audit-sha256", "unused",
        "--binding", "absent", "--binding-sha256", "unused", "--prior-binding", "absent",
        "--prior-binding-sha256", "unused", "--output", str(output)])
    with pytest.raises(ValueError, match="Output already exists"):
        reader.main()
    assert output.read_text() == "retain me"
