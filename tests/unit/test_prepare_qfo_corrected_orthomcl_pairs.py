import io
from itertools import combinations
import random
from pathlib import Path

import pytest

from benchmark_tools import prepare_qfo_corrected_orthomcl_pairs as module


def admission():
    return {"status": "corrected_orthomcl_native_outputs_admitted", "accuracy_admitted": False,
            "publication_ready": False, "pair_semantics": module.SEMANTICS,
            "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "180", "ReqMem": "900G"},
            "content": {"input_proteins": 984137, "input_species": 78, "cross_species_clique_pairs": 30,
                        "final_groups": 12, "grouped_proteins": 41, "ungrouped_input_proteins": 984096}}


def test_admission():
    value = admission()
    assert module.validate_admission(value) == value["content"]


@pytest.mark.parametrize("key,value", [("status", "failed"), ("accuracy_admitted", True),
    ("publication_ready", True), ("pair_semantics", "pre-clustering graph edges")])
def test_wrong_admission(key, value):
    data = admission()
    data[key] = value
    with pytest.raises(ValueError):
        module.validate_admission(data)


@pytest.mark.parametrize("key,value", [("input_proteins", 976504), ("input_species", 77),
    ("cross_species_clique_pairs", -1), ("cross_species_clique_pairs", True),
    ("final_groups", -1), ("grouped_proteins", -1), ("ungrouped_input_proteins", 0)])
def test_wrong_counts(key, value):
    data = admission()
    data["content"][key] = value
    with pytest.raises(ValueError):
        module.validate_admission(data)


@pytest.mark.parametrize("key", ["State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem"])
def test_wrong_scheduler(key):
    data = admission()
    data["scheduler"][key] = "wrong"
    with pytest.raises(ValueError):
        module.validate_admission(data)


def fixture(tmp_path):
    groups = tmp_path / "groups"
    groups.write_text("ORTHOMCL0(3 genes,2 taxa): sp|B|b(x) sp|A|a(x) sp|C|c(y)\n"
                      "ORTHOMCL2(2 genes,1 taxa): D(y) E(y)\n")
    owners = {"sp|B|b": "x", "sp|A|a": "x", "sp|C|c": "y", "D": "y", "E": "y", "unassigned": "x"}
    return groups, owners


def test_cross_species_only_and_unassigned_coverage(tmp_path):
    groups, owners = fixture(tmp_path)
    stream = io.StringIO()
    assert module.write_pairs(groups, owners, stream) == {
        "total_pairs": 2, "final_groups": 2, "grouped_proteins": 5, "ungrouped_input_proteins": 1}
    assert stream.getvalue() == "A\tC\nB\tC\n"


@pytest.mark.parametrize("old,new", [
    ("sp|B|b(x)", "sp|A|a(x)"), ("sp|B|b(x)", "unknown(x)"),
    ("D(y)", "sp|B|b(x)"), ("ORTHOMCL2", "ORTHOMCL0"),
    ("D(y) E(y)", "D(y)"), ("D(y)", "bad"),
])
def test_invalid_partition(tmp_path, old, new):
    groups, owners = fixture(tmp_path)
    groups.write_text(groups.read_text().replace(old, new))
    with pytest.raises(ValueError):
        module.write_pairs(groups, owners, io.StringIO())


def test_collision_even_in_ungrouped_inputs(tmp_path):
    groups, owners = fixture(tmp_path)
    owners["tr|B|other"] = "z"
    with pytest.raises(ValueError, match="injective"):
        module.write_pairs(groups, owners, io.StringIO())


def test_empty_predictions_do_not_invent_pairs(tmp_path):
    groups, owners = fixture(tmp_path)
    groups.write_text("")
    stream = io.StringIO()
    counts = module.write_pairs(groups, owners, stream)
    assert counts == {"total_pairs": 0, "final_groups": 0, "grouped_proteins": 0, "ungrouped_input_proteins": 6}
    assert stream.getvalue() == ""


@pytest.mark.parametrize("seed", range(10))
def test_random_disjoint_groups_against_brute_force(tmp_path, seed):
    rng, owners, expected, lines = random.Random(seed), {}, set(), []
    for cid in range(30):
        genes = [f"G{cid}_{i}" for i in range(rng.randint(2, 16))]
        for gene in genes:
            owners[gene] = f"S{rng.randrange(5)}"
        expected.update(tuple(sorted((a, b))) for a, b in combinations(genes, 2) if owners[a] != owners[b])
        lines.append(f"ORTHOMCL{cid}({len(genes)} genes,{len({owners[g] for g in genes})} taxa): "
                     + " ".join(f"{g}({owners[g]})" for g in genes))
    owners["unassigned"] = "S0"
    path = tmp_path / "groups"
    path.write_text("\n".join(lines) + "\n")
    output = io.StringIO()
    counts = module.write_pairs(path, owners, output)
    pairs = [tuple(line.split("\t")) for line in output.getvalue().splitlines()]
    assert set(pairs) == expected
    assert len(pairs) == len(set(pairs)) == counts["total_pairs"]
    assert counts["ungrouped_input_proteins"] == 1


def test_real_retained_native_fixture():
    root = Path(__file__).resolve().parents[2] / "benchmarks/work/orthomcl_staged_native_inference_probe_v2_20260918"
    if not root.exists():
        pytest.skip("Retained native fixture not available")
    owners = module.load_species(root / "inputs/all.gg")
    groups = root / "tool/Sep_18/all_orthomcl.out"
    stream = io.StringIO()
    counts = module.write_pairs(groups, owners, stream)
    expected = {tuple(sorted(pair)) for _, genes in module.iter_orthomcl(groups)
                for pair in combinations(genes, 2) if owners[pair[0]] != owners[pair[1]]}
    observed = [tuple(line.split("\t")) for line in stream.getvalue().splitlines()]
    assert set(observed) == expected
    assert counts == {"total_pairs": 30, "final_groups": 12, "grouped_proteins": 41, "ungrouped_input_proteins": 1}
    assert len(observed) == len(set(observed)) == 30


def test_prepare_requires_scheduled_allocation(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, tmp_path / "admission", "sha", 21752)


def test_prepare_rejects_pending_admission_before_output(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "65536")
    monkeypatch.setattr(module.os, "uname", lambda: type("Host", (), {"nodename": "bizon"})())
    monkeypatch.setattr(module, "verify_runtime", lambda root: {})
    def pending(*args):
        raise ValueError("Pending admission")
    monkeypatch.setattr(module, "completed", pending)
    with pytest.raises(ValueError, match="Pending"):
        module.prepare(tmp_path, tmp_path / "admission", "sha", 21752)
    assert not (tmp_path / "benchmarks").exists()
