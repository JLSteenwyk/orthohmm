import copy
import json
from pathlib import Path

import pytest

import benchmark_tools.run_corrected_swiss_strata as module
from benchmark_tools.bootstrap_corrected_swiss_strata import METHODS

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def fixture():
    # Real input-only features; all prediction counts below are synthetic.
    strata = json.loads((RESULTS / "corrected_swiss_sequence_strata_20260918.json").read_text())
    rows = [dict(family=f, represented_genes=genes,
                 counts_without_prior=dict(TP=2, FN=1, FP=1, TN=2))
            for f, genes in strata["family_memberships"].items()]
    shared = dict(reference=dict(sha256=module.REFERENCE_SHA), shared_represented_genes={},
                  reference_relation_count=10765)
    factorial = dict(**shared, status="corrected_qfo_factorial_swiss_counts_verified", checked_inputs=[],
                     cells=[dict(cell=f"p{p}_c{c}_r{r}", families=copy.deepcopy(rows))
                            for p in (0, 1) for c in (0, 1) for r in (0, 1)])
    comparator = dict(**shared, status="corrected_comparator_swiss_counts_verified", checked_records=[],
                      method="orthofinder_full", method_key="orthofinder_3_1_5_full", families=copy.deepcopy(rows))
    return factorial, comparator, strata


def test_real_frozen_input_bins_and_exact_cell_selection():
    factorial, comparator, strata = fixture()
    counts, membership = module.assemble(factorial, comparator, strata)
    assert set(counts) == set(METHODS)
    assert list(membership) == sorted(strata["family_memberships"])
    assert {f for f, b in membership.items() if b == "lower"} == set(strata["primary_strata"]["lower_entropy"])
    assert list(counts[METHODS[0]].values()) == factorial["cells"][4]["families"]
    assert list(counts[METHODS[1]].values()) == factorial["cells"][7]["families"]


def test_cell_labels_match_upstream_count_auditor_contract():
    from benchmark_tools.audit_qfo_factorial_swiss import CELLS as upstream_cells
    factorial, comparator, strata = fixture()
    assert [r["cell"] for r in factorial["cells"]] == list(upstream_cells)
    assert upstream_cells[4] == "p1_c0_r0" and upstream_cells[7] == "p1_c1_r1"
    module.assemble(factorial, comparator, strata)
    for row in factorial["cells"]:
        row["cell"] = row["cell"].replace("_", "")
    with pytest.raises(ValueError, match="inventory/order"):
        module.assemble(factorial, comparator, strata)


@pytest.mark.parametrize("change", ["historical", "sequence_only", "reference", "overlap", "relations",
                                    "cell_order", "duplicate_cell", "membership", "bin", "cutoff", "scored"])
def test_mixed_or_changed_evidence_rejected(change):
    f, c, s = fixture()
    if change == "historical":
        f["status"] = "qfo_factorial_swiss_counts_verified"
    elif change == "sequence_only":
        c["method"] = "orthofinder_sequence_only"
    elif change == "reference":
        c["reference"] = dict(sha256="wrong")
    elif change == "overlap":
        f["shared_represented_genes"] = {"g": ["a", "b"]}
    elif change == "relations":
        f["reference_relation_count"] -= 1
    elif change == "cell_order":
        f["cells"].reverse()
    elif change == "duplicate_cell":
        f["cells"][4] = f["cells"][0]
    elif change == "membership":
        c["families"][0]["represented_genes"][0] = "foreign"
    elif change == "bin":
        s["primary_strata"]["lower_entropy"].pop()
    elif change == "cutoff":
        s["median_family_entropy_cutoff"] += .01
    else:
        s["prediction_statistics_evaluated"] = True
    with pytest.raises(ValueError):
        module.assemble(f, c, s)


def setup_run(tmp_path, monkeypatch):
    factorial, comparator, strata = fixture()
    evidence = tmp_path / "evidence.json"
    evidence.write_text("{}")
    identity = module.record(evidence)
    for key in ("source", "helper", "descriptor_update_protocol", "native_sequence_audit"):
        strata[key] = identity
    strata["inputs"] = [identity]
    source = tmp_path / "strata.json"
    source.write_text(json.dumps(strata))
    monkeypatch.setattr(module, "STRATA_SHA", module.record(source)["sha256"])
    calls = []
    def audit_f(*args):
        calls.append(("factorial", args))
        return factorial
    def audit_c(*args):
        calls.append(("comparator", args))
        return comparator
    monkeypatch.setattr(module, "audit_factorial", audit_f)
    monkeypatch.setattr(module, "audit_comparator", audit_c)
    real_bootstrap = module.bootstrap
    monkeypatch.setattr(module, "bootstrap", lambda c, m: real_bootstrap(c, m, replicates=100))
    output = tmp_path / "result.json"
    args = (evidence, identity["sha256"], evidence, identity["sha256"], evidence, source,
            RESULTS / "CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md", output)
    return args, calls


def test_driver_reconstructs_counts_and_retains_provenance(tmp_path, monkeypatch):
    args, calls = setup_run(tmp_path, monkeypatch)
    result = module.run(*args)
    assert [c[0] for c in calls] == ["factorial", "comparator"]
    assert calls[0][1] == (args[0], args[1], args[4])
    assert calls[1][1] == (args[2], args[3], args[4])
    assert result["scientific_inputs_admitted"] is True
    assert result["publication_ready"] is False
    assert result["method_bindings"]["high_sensitivity"] == "corrected replay p1c0r0"
    assert json.loads(args[-1].read_text()) == result
    with pytest.raises(FileExistsError):
        module.run(*args)


@pytest.mark.parametrize("change", ["helper", "protocol", "strata", "postcheck", "audit"])
def test_driver_fails_closed_without_result(tmp_path, monkeypatch, change):
    args, calls = setup_run(tmp_path, monkeypatch)
    if change == "helper":
        monkeypatch.setattr(module, "SOURCES", {"bootstrap_corrected_swiss_strata.py": "wrong"})
    elif change == "protocol":
        monkeypatch.setattr(module, "PROTOCOL_SHA", "wrong")
    elif change == "strata":
        monkeypatch.setattr(module, "STRATA_SHA", "wrong")
    elif change == "audit":
        def fail(*unused):
            raise ValueError("raw evidence failed")
        monkeypatch.setattr(module, "audit_factorial", fail)
    else:
        real = module.bootstrap
        def mutate(c, m):
            result = real(c, m)
            args[0].write_text("changed during computation")
            return result
        monkeypatch.setattr(module, "bootstrap", mutate)
    with pytest.raises((ValueError, RuntimeError)):
        module.run(*args)
    assert not args[-1].exists()
