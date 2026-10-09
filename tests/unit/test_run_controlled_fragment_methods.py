from copy import deepcopy
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from benchmark_tools.controlled_fragment_observations import CONDITION, SEEDS, INFERENCE_METHODS, observe, verify_observation
from benchmark_tools.prepare_controlled_fragment_observations import record
from benchmark_tools.run_controlled_fragment_methods import allocation, observed_inputs, run


@pytest.fixture
def dataset(tmp_path):
    parent = tmp_path / "parent/input"
    observed = tmp_path / "observed/input"
    parent.mkdir(parents=True)
    observed.mkdir(parents=True)
    sequences = {f"g{i}": "ACDEFGHIKLMNPQRSTVWY" for i in range(10)}
    owners = {g: "A" if i % 2 else "B" for i, g in enumerate(sequences)}
    for species in ("A", "B"):
        SeqIO.write([SeqRecord(Seq(seq), id=g, description="") for g, seq in sequences.items() if owners[g] == species], parent / (species + ".fasta"), "fasta")
    truth = {"extant_genes": 10, "species": ["A", "B"], "families": {"one": sorted(sequences)},
             "ortholog_pairs": [["g0", "g1"]], "ortholog_pair_count": 1,
             "prepared_inputs": [{"path": "input/" + p.name, "bytes": p.stat().st_size, "sha256": record(p)["sha256"]} for p in sorted(parent.iterdir())]}
    old_truth = parent.parent / "truth.json"
    new_truth = observed.parent / "parent_truth.json"
    old_truth.write_text(json.dumps(truth))
    new_truth.write_bytes(old_truth.read_bytes())
    transformed, coordinates = observe(sequences, owners, SEEDS[0])
    for species in ("A", "B"):
        SeqIO.write([SeqRecord(Seq(seq), id=g, description="") for g, seq in transformed.items() if owners[g] == species], observed / (species + ".fasta"), "fasta")
    coordinate_path = observed.parent / "coordinates.json"
    coordinate_path.write_text(json.dumps(coordinates))
    return {"input": str(observed), "truth": str(new_truth), "seed": SEEDS[0],
            "verified_inputs": {"status": "ready", "truth": record(new_truth), "inputs": [record(p) for p in sorted(observed.iterdir())]},
            "parent_inputs": {"status": "ready", "truth": record(old_truth), "inputs": [record(p) for p in sorted(parent.iterdir())]},
            "coordinates": record(coordinate_path), "verification": verify_observation(sequences, transformed, owners, coordinates, SEEDS[0])}


def test_readback_validates_real_written_fasta_and_unchanged_truth(dataset):
    assert observed_inputs(dataset) == dataset["verified_inputs"]


@pytest.mark.parametrize("fault", ["extra_file", "sequence", "truth", "coordinates", "verification", "readiness", "truth_path"])
def test_disk_readback_rejects_mutations_even_when_hash_is_rebound(dataset, fault, tmp_path):
    if fault == "extra_file":
        (Path(dataset["input"]) / "extra.txt").write_text("extra")
    elif fault == "sequence":
        item = dataset["verified_inputs"]["inputs"][0]
        path = Path(item["absolute_path"])
        path.write_text(path.read_text().replace("ACD", "XXX"))
        dataset["verified_inputs"]["inputs"][0] = record(path)
    elif fault == "truth":
        path = Path(dataset["truth"])
        path.write_text(path.read_text() + "\n")
        dataset["verified_inputs"]["truth"] = record(path)
    elif fault == "coordinates":
        path = Path(dataset["coordinates"]["absolute_path"])
        rows = json.loads(path.read_text())
        rows[0]["start"] += 1
        path.write_text(json.dumps(rows))
        dataset["coordinates"] = record(path)
    elif fault == "verification":
        dataset["verification"]["fragment_genes"] = 0
    elif fault == "readiness":
        dataset["verified_inputs"]["status"] = "other"
    else:
        dataset["truth"] = str(tmp_path / "other")
    with pytest.raises(ValueError):
        observed_inputs(dataset)


@pytest.mark.parametrize("cpus,memory,job", [(0, 16384, "1"), (4, 16384, "1"), (16, 1024, "1"), (16, 16384, "")])
def test_scheduler_limits_required(monkeypatch, cpus, memory, job):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", str(cpus))
    monkeypatch.setenv("SLURM_MEM_PER_NODE", str(memory))
    monkeypatch.setenv("SLURM_JOB_ID", job)
    with pytest.raises(ValueError, match="scheduler allocation"):
        allocation()


def test_scheduler_records_shared_allocation_not_isolation(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "16384")
    monkeypatch.setenv("SLURM_JOB_ID", "1")
    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.subprocess.check_output", lambda *a, **k: "actual allocation")
    result = allocation()
    assert result["scontrol"] == "actual allocation" and result["cpus_per_task"] == 16
    assert result["scope"].endswith("not node isolation")


@pytest.mark.parametrize("fault", [None, "native_failure", "preflight_failure", "interrupt"])
def test_panel_runs_in_exact_order_without_retrying_any_failure(tmp_path, monkeypatch, fault):
    panel = {"datasets": [{"label": f"{CONDITION}_{s}", "seed": s} for s in SEEDS]}
    manifest, runtime_path = tmp_path / "manifest.json", tmp_path / "runtime.json"
    manifest.write_text("{}")
    runtime_path.write_text("{}")
    called = []
    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.load_panel", lambda *a: (panel, {}, {}))
    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.allocation", lambda: {"job_id": "1"})
    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.subprocess.check_output", lambda *a, **k: "frozen")
    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.read_frozen", lambda *a: {})
    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.verify_environment", lambda *a: None)

    def inputs(d):
        if fault == "preflight_failure" and d["seed"] == SEEDS[2]:
            raise ValueError("affected dataset only")
        return {"status": "ready"}

    def execute(d, order, env, evidence, verified, provenance, runtime):
        called.append((d["seed"], tuple(order)))
        if fault == "interrupt" and d["seed"] == SEEDS[2]:
            raise KeyboardInterrupt("interrupted")
        evidence.mkdir(parents=True)
        (evidence / "status.json").write_text("{}")
        return {"failed_methods": [INFERENCE_METHODS[0]] if fault == "native_failure" and d["seed"] == SEEDS[2] else []}

    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.observed_inputs", inputs)
    monkeypatch.setattr("benchmark_tools.run_controlled_fragment_methods.execute", execute)
    destination = tmp_path / "inference"
    args = (tmp_path, manifest, "sha", runtime_path, "runtime_sha", destination)
    if fault == "interrupt":
        with pytest.raises(KeyboardInterrupt):
            run(*args)
        result = json.loads((destination / "panel_status.json").read_text())
        assert result["status"] == "interrupted_or_integrity_failure"
        assert [s for s, _ in called] == list(SEEDS[:3])
    else:
        result = run(*args)
        expected = [s for s in SEEDS if not (fault == "preflight_failure" and s == SEEDS[2])]
        assert [s for s, _ in called] == expected
        assert len(result["datasets"]) == 10 and result["accuracy_evaluated"] is False
        assert result["status"] == "finished_pending_native_validation"
        if fault == "preflight_failure":
            assert result["datasets"][2]["status"] == "failed"
        elif fault == "native_failure":
            assert result["datasets"][2]["failed_methods"] == [INFERENCE_METHODS[0]]
    assert all(order == INFERENCE_METHODS for _, order in called)
    assert len({s for s, _ in called}) == len(called)
    with pytest.raises(FileExistsError, match="no automatic restart"):
        run(*args)
