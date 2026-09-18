from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import run_qfo_sequence_uncertainty as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_bootstrap_qfo_sequence import fixture


def setup(tmp_path, monkeypatch):
    directory = Path(module.__file__).parent
    protocol = directory / "results/QFO_SEQUENCE_UNCERTAINTY_PROTOCOL_20260918.md"
    anchor = tmp_path / "anchor.json"
    anchor.write_text("{}\n")
    data = fixture()
    data.update(source=record(directory / "audit_qfo_sequence_swiss.py"),
        admission_inventory=record(anchor), baseline_audit=record(anchor),
        checked_inputs=[record(anchor)], helpers=[record(directory / name) for name in module.SOURCES])
    counts = tmp_path / "counts.json"
    counts.write_text(json.dumps(data))
    calls = []

    def audit(path, digest, baseline):
        calls.append((path, digest, baseline))
        return deepcopy(data)

    monkeypatch.setattr(module, "audit", audit)
    return data, counts, protocol, anchor, calls


def test_fixed_production_controls_and_source_reconstruction(tmp_path, monkeypatch):
    data, counts, protocol, anchor, calls = setup(tmp_path, monkeypatch)
    # Run the actual 100,000-draw production defaults on synthetic counts.
    result = module.run(counts, record(counts)["sha256"], protocol, tmp_path / "out.json", tmp_path / "out.md")
    assert calls == [(anchor, record(anchor)["sha256"], anchor)]
    assert result["replicates"] == 100000 and result["seed"] == 20260923
    assert result["uncertainty_admitted"] is True and result["publication_ready"] is False
    assert result["multiplicity_endpoints"] == 6
    assert len(result["comparisons"]) == 2
    assert json.loads((tmp_path / "out.json").read_text())["comparisons"] == result["comparisons"]
    rendered = (tmp_path / "out.md").read_text()
    assert rendered.count("| all_hits - p0_c0_r0 |") == 3
    assert rendered.count("| top100 - p0_c0_r0 |") == 3
    assert "Development-exposed" in rendered
    with pytest.raises(FileExistsError):
        module.run(counts, record(counts)["sha256"], protocol, tmp_path / "out.json", tmp_path / "new.md")


@pytest.mark.parametrize("problem", ["protocol", "source", "hash", "raw", "reconstruction", "helper"])
def test_changed_evidence_rejected_without_output(tmp_path, monkeypatch, problem):
    data, counts, protocol, anchor, calls = setup(tmp_path, monkeypatch)
    digest = record(counts)["sha256"]
    if problem == "protocol":
        protocol = anchor
    elif problem == "source":
        data["source"] = record(anchor)
        counts.write_text(json.dumps(data))
        digest = record(counts)["sha256"]
    elif problem == "hash":
        digest = "0" * 64
    elif problem == "raw":
        anchor.write_text("changed\n")
    elif problem == "reconstruction":
        data["variants"][0]["families"][0]["counts_without_prior"]["TP"] += 1
    else:
        monkeypatch.setattr(module, "SOURCES", {**module.SOURCES, "bootstrap_qfo_sequence.py": "0" * 64})
    with pytest.raises(ValueError):
        module.run(counts, digest, protocol, tmp_path / "out.json", tmp_path / "out.md")
    assert not (tmp_path / "out.json").exists()
    assert not (tmp_path / "out.md").exists()


@pytest.mark.parametrize("kind", ["same", "symlink"])
def test_output_collisions(tmp_path, monkeypatch, kind):
    _, counts, protocol, _, _ = setup(tmp_path, monkeypatch)
    output, markdown = tmp_path / "out.json", tmp_path / "out.md"
    if kind == "same":
        markdown = output
    else:
        output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        module.run(counts, record(counts)["sha256"], protocol, output, markdown)


def test_frozen_hashes():
    directory = Path(module.__file__).parent
    assert record(directory / "results/QFO_SEQUENCE_UNCERTAINTY_PROTOCOL_20260918.md")["sha256"] == module.PROTOCOL_SHA
    for name, digest in module.SOURCES.items():
        assert record(directory / name)["sha256"] == digest
