from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import bootstrap_qfo_corrected_factorial as module
from benchmark_tools.bootstrap_qfo_factorial import bootstrap
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_bootstrap_qfo_factorial import fixture


def corrected():
    return {**fixture(), "status": "corrected_qfo_factorial_swiss_counts_verified",
            "uncertainty_admitted": False, "publication_ready": False}


def test_numerically_identical_protocol_without_input_mutation():
    data = corrected()
    saved = deepcopy(data)
    result = module.calculate(data, replicates=100, seed=57)
    original = bootstrap(fixture(), replicates=100, seed=57)
    assert data == saved
    assert result["status"] == "paired_corrected_qfo_factorial_swiss_intervals"
    for key in ("point_estimates", "comparisons", "seed", "replicates", "multiplicity_endpoints", "rng"):
        assert result[key] == original[key]


@pytest.mark.parametrize("mutation", ["historical", "admitted", "publication", "missing", "overlap", "truth"])
def test_invalid_input(mutation):
    data = corrected()
    if mutation == "historical":
        data["status"] = "qfo_factorial_swiss_counts_verified"
    elif mutation == "admitted":
        data["uncertainty_admitted"] = True
    elif mutation == "publication":
        data["publication_ready"] = True
    elif mutation == "missing":
        data["cells"].pop()
    elif mutation == "overlap":
        data["shared_represented_genes"] = {"gene": ["a", "b"]}
    else:
        data["cells"][1]["families"][0]["counts_without_prior"]["TP"] += 1
    with pytest.raises(ValueError):
        module.calculate(data, replicates=100, seed=57)


def test_frozen_implementation_hashes():
    directory = Path(module.__file__).parent
    assert record(directory / "bootstrap_qfo_factorial.py")["sha256"] == module.ENGINE_SHA
    assert record(directory / "audit_qfo_corrected_factorial_swiss.py")["sha256"] == module.AUDITOR_SHA
    assert record(directory / "results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md")["sha256"] == module.CORRECTED_PROTOCOL_SHA


def test_file_bound_run(tmp_path, monkeypatch):
    directory = Path(module.__file__).parent
    protocol = directory / "results/QFO_FACTORIAL_PROTOCOL_20260917.md"
    corrected_protocol = directory / "results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md"
    anchor = tmp_path / "anchor.json"
    anchor.write_text("{}\n")
    data = corrected()
    data.update(source=record(directory / "audit_qfo_corrected_factorial_swiss.py"),
                admission_inventory=record(anchor), baseline_audit=record(anchor),
                checked_inputs=[record(anchor)], helpers=[])
    counts = tmp_path / "counts.json"
    counts.write_text(json.dumps(data))
    digest = record(counts)["sha256"]
    original = module.calculate
    calls = []
    def small_run(report, **kwargs):
        calls.append(kwargs)
        return original(report, replicates=100, seed=20260922)
    monkeypatch.setattr(module, "calculate", small_run)
    result = module.run(counts, digest, protocol, corrected_protocol, tmp_path / "result.json", tmp_path / "result.md")
    assert calls == [{}]  # Production uses fixed defaults, not caller-supplied controls.
    assert "Corrected-Release" in (tmp_path / "result.md").read_text()
    assert result["counts"]["sha256"] == digest
    with pytest.raises(FileExistsError):
        module.run(counts, digest, protocol, corrected_protocol, tmp_path / "result.json", tmp_path / "another.md")
    anchor.write_text("changed")
    with pytest.raises(ValueError):
        module.run(counts, digest, protocol, corrected_protocol, tmp_path / "bad.json", tmp_path / "bad.md")
    assert not (tmp_path / "bad.json").exists()
