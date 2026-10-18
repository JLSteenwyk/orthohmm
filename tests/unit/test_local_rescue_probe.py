import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from benchmark_tools import run_local_rescue_probe as driver
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import PRODUCTS
from benchmark_tools.verify_frozen_phylogeny_install import fixture
from orthohmm.accuracy import write_accuracy_checkpoint


def case_fixture(root, dataset="OrthoBench", active=True):
    root.mkdir(parents=True)
    fasta, phylo = root / "input", root / "baseline/orthohmm_phylogeny"
    fasta.mkdir()
    phylo.mkdir(parents=True)
    sequences = fixture()
    for taxon, genes in sequences.items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    names = sorted(gene for genes in sequences.values() for gene in genes)
    species = np.array([int(gene[1]) - 1 for gene in names], dtype=np.int32)
    groups = [[gene for gene in names if gene.endswith(f"family{index}")] for index in range(3)]
    groups.append([gene for gene in names if gene.endswith("duplicate")])
    working = root / "candidate"
    working.mkdir()
    candidate, constraints = working / "original.txt", working / "constraints.json"
    candidate.write_text("".join(" ".join(group) + "\n" for group in groups))
    constraints.write_text("[]\n")
    links = [("S1_family2", "S2_family2_duplicate"), ("S2_family2", "S3_family2_duplicate")]
    links += [(b, a) for a, b in links]
    if not active:
        links = []
    index = {gene: number for number, gene in enumerate(names)}
    checkpoint = write_accuracy_checkpoint(str(root / "initial"), names, species,
        [index[a] for a, b in links], [index[b] for a, b in links], [9.] * len(links))
    rows = [f"{a}\t{a[:2]}\t{b}\t{b[:2]}\thigh\n" for group in groups
            for i, a in enumerate(group) for b in group[i + 1:]]
    (phylo / PRODUCTS[2]).write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\n" + "".join(rows))
    for name in PRODUCTS[:2]:
        (phylo / name).write_text("unused_in_label_free_fixture\n")
    (phylo / "species_tree.rooted.nwk").write_text("[&R] ((S1,S2),(S3,S4));\n")
    return {"dataset": dataset, "universe": [16, 4], "fasta_directory": str(fasta),
            "fastas": [record(path) for path in sorted(fasta.iterdir())], "checkpoint": str(checkpoint),
            "checkpoint_manifest": record(checkpoint / "manifest.json"), "candidate": record(candidate),
            "constraints": record(constraints), "baseline_output": str(phylo.parent),
            "golden": {path.name: record(path) for path in sorted(phylo.iterdir())},
            "aligner": shutil.which("mafft") or "mafft", "tree_builder": shutil.which("FastTree") or "FastTree"}


def execution_plan(tmp_path):
    cases = [case_fixture(tmp_path / dataset, dataset) for dataset in driver.UNIVERSES]
    inputs = []
    for case in cases:
        case["universe"] = list(driver.UNIVERSES[case["dataset"]])
        inputs.extend([*case["fastas"], case["candidate"], case["constraints"], case["checkpoint_manifest"], *case["golden"].values()])
    output = tmp_path / "output"
    output.mkdir()
    plan = {"status": "prepared_label_free_probe_unrun", "source_root": str(Path(driver.__file__).resolve().parents[1]),
            "source_revision": "fixture", "output": str(output), "cpu": 16, "cases": cases,
            "inputs": inputs, "sources": [], "runtime": {"files": []}, "tool_path": "/usr/bin",
            "scope": "test coordinator only"}
    path = output / "execution_plan.json"
    path.write_text(json.dumps(plan))
    return plan, path


def test_fixed_parameters_still_match_the_prespecified_protocol():
    protocol = Path(driver.__file__).resolve().parents[1] / driver.PROTOCOL
    assert json.loads(protocol.read_text())["parameters"] == driver.PARAMETERS


@pytest.mark.parametrize("changed", [None, "species_tree.rooted.nwk", PRODUCTS[1]])
def test_original_qfo_scoring_failure_requires_independently_admitted_unchanged_native_products(tmp_path, monkeypatch, changed):
    repo = tmp_path / "repo"
    canonical = repo / driver.CANONICAL
    configurations, rows, admitted = [], [], {}
    for dataset, universe in driver.UNIVERSES.items():
        case = case_fixture(tmp_path / dataset, dataset)
        native = canonical / dataset / "canonical/native"
        shutil.copytree(Path(case["baseline_output"]), native)
        golden = {path.name: record(path) for path in (native / "orthohmm_phylogeny").iterdir()}
        admitted[dataset] = {"coverage": {"genes": universe[0], "species": universe[1]},
                             "products": {name: golden[name] for name in PRODUCTS},
                             "checked": [golden["species_tree.rooted.nwk"]]}
        configurations.append({key: case[key] for key in ("dataset", "fastas", "fasta_directory", "aligner", "tree_builder")})
        rows.append({"dataset": dataset, "label": "canonical", "status": "failed" if dataset == "QfO" else "verified",
            "ordering": {"manifest": case["checkpoint_manifest"]},
            "candidate": {"candidate_partition": case["candidate"], "membership_constraints": case["constraints"]}})
    path = canonical / "execution_plan.json"
    path.write_text(json.dumps({"cases": configurations}))
    (canonical / "results.json").write_text(json.dumps({"execution_plan": record(path), "arms": rows}))
    evidence = repo / driver.ADMISSION
    evidence.parent.mkdir(parents=True)
    evidence.write_text("{}")
    monkeypatch.setattr(driver, "verified_cell", lambda summary, key, dataset: admitted[dataset])
    if changed:
        altered = canonical / "QfO/canonical/native/orthohmm_phylogeny" / changed
        altered.write_text(altered.read_text() + "\n")
        with pytest.raises(ValueError, match="admission"):
            driver.selected_cases(repo)
    else:
        cases, pins = driver.selected_cases(repo)
        assert [case["dataset"] for case in cases] == list(driver.UNIVERSES)
        assert cases[1]["original_producer_status"] == "failed" and len(pins) == 3


@pytest.mark.parametrize("kind", ["order", "budget", "universe", "pin", "directory", "checkpoint", "golden"])
def test_execution_plan_rejects_unpinned_paths_or_changed_contrast(tmp_path, kind):
    plan, _ = execution_plan(tmp_path)
    driver.validate_plan(plan)
    if kind == "order":
        plan["cases"].reverse()
    elif kind == "budget":
        plan["cpu"] = 32
    elif kind == "universe":
        plan["cases"][0]["universe"][0] -= 1
    elif kind == "pin":
        plan["inputs"].pop()
    elif kind == "directory":
        plan["cases"][0]["fasta_directory"] = str(tmp_path / "other")
    elif kind == "checkpoint":
        plan["cases"][0]["checkpoint"] = str(tmp_path / "other")
    else:
        plan["cases"][0]["baseline_output"] = str(tmp_path / "other")
    with pytest.raises(ValueError):
        driver.validate_plan(plan)


@pytest.mark.parametrize("statuses,expected,visited,admitted", [
    (["analyzable_contexts_found"], "applicability_found_accuracy_evaluation_pending", 1, True),
    (["analyzable_contexts_found"], "analyzable_but_zero_admitted_no_accuracy_claim", 1, False),
    (["no_analyzable_contexts", "analyzable_contexts_found"], "applicability_found_accuracy_evaluation_pending", 2, True),
    (["no_analyzable_contexts"] * 3, "all_three_datasets_inactive_no_accuracy_claim", 3, False),
    (["failed"], "failed_probe_no_fallback_on_failure", 1, False),
])
def test_coordinator_priority_fallback_and_failure_are_not_interchangeable(tmp_path, monkeypatch, statuses, expected, visited, admitted):
    plan, path = execution_plan(tmp_path)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(driver, "runtime", lambda: plan["runtime"])
    monkeypatch.setattr(driver, "host_activity", lambda: {"fixture": True})
    calls = []
    def stage(command, output, name, source, env):
        dataset = command[-1]
        coordinator = json.loads((output / "results.json").read_text())
        assert coordinator["status"] == "running" and coordinator["active_dataset"] == dataset
        assert env["NUMBA_NUM_THREADS"] == "16" and env["CUDA_VISIBLE_DEVICES"] == ""
        folder = output / dataset
        folder.mkdir()
        status = statuses[len(calls)]
        calls.append(dataset)
        found = status == "analyzable_contexts_found"
        (folder / "probe.json").write_text(json.dumps({"status": status, "tree_contexts": int(found), "accepted_contexts": int(found and admitted)}))
        return {"returncode": int(status == "failed"), "command": command}
    monkeypatch.setattr(driver, "timed_stage", stage)
    result = driver.run(path)
    assert result["status"] == expected
    assert calls == list(driver.UNIVERSES)[:visited]
    assert result["unrun_datasets"] == list(driver.UNIVERSES)[visited:]
    assert result["active_dataset"] is None and not result["accuracy_evaluated"]
    with pytest.raises(ValueError, match="preserve"):
        driver.run(path)


def test_context_pair_reader_keeps_only_high_pairs_and_rejects_partial_context_ownership(tmp_path):
    path = tmp_path / "confidence.tsv"
    header = "gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\n"
    contexts = [{"context_id": "LocalContext000001", "genes_a": ["a1", "a2"], "genes_b": ["b1", "b2"]}]
    species = {"a1": "S1", "a2": "S2", "b1": "S1", "b2": "S2", "x1": "S1", "x2": "S2"}
    path.write_text(header + "a1\tS1\ta2\tS2\thigh\nb1\tS1\tb2\tS2\tmedium\nx1\tS1\tx2\tS2\thigh\n")
    assert driver.context_baseline_pairs(path, contexts, species) == {"LocalContext000001": [("a1", "a2")]}
    path.write_text(header + "a1\tS1\tx2\tS2\thigh\n")
    with pytest.raises(ValueError, match="ownership"):
        driver.context_baseline_pairs(path, contexts, species)
    with pytest.raises(ValueError, match="overlap"):
        driver.context_baseline_pairs(path, contexts + contexts, species)


def test_changed_input_during_probe_is_a_terminal_verification_failure(tmp_path, monkeypatch):
    plan, path = execution_plan(tmp_path)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(driver, "runtime", lambda: plan["runtime"])
    monkeypatch.setattr(driver, "host_activity", lambda: {})
    def changed(command, output, name, source, env):
        folder = output / command[-1]
        folder.mkdir()
        (folder / "probe.json").write_text(json.dumps({"status": "analyzable_contexts_found", "tree_contexts": 1, "accepted_contexts": 1}))
        Path(plan["cases"][0]["candidate"]["path"]).write_text("changed\n")
        return {"returncode": 0}
    monkeypatch.setattr(driver, "timed_stage", changed)
    result = driver.run(path)
    assert result["status"] == "failed_final_verification" and result["error_type"] == "ValueError"
    assert json.loads((Path(plan["output"]) / "results.json").read_text()) == result


def test_inactive_probe_preserves_exact_candidates_and_skips_all_local_scoring(tmp_path, monkeypatch):
    case = case_fixture(tmp_path / "case", active=False)
    def forbidden(*args, **kwargs):
        raise AssertionError("No nominated context must not score profiles or infer trees")
    monkeypatch.setattr(driver, "load_global_sequence_database", forbidden)
    result = driver.probe_case({"cpu": 2}, case, tmp_path / "probe")
    assert result["status"] == "no_analyzable_contexts"
    assert result["tree_contexts"] == result["accepted_contexts"] == 0
    assert not result["accuracy_evaluated"] and not result["final_native_inference"]
    assert Path(result["candidate"]["path"]).read_bytes() == Path(case["candidate"]["path"]).read_bytes()
    assert result["effective_partition"]["original_bytes_preserved"]


@pytest.mark.skipif(not shutil.which("mafft") or not shutil.which("FastTree"), reason="Requires real tree tools")
def test_real_probe_connects_nomination_profiles_fresh_tree_and_effective_partition(tmp_path):
    case = case_fixture(tmp_path / "case")
    before = {pin["path"]: record(pin["path"]) for pin in [case["candidate"], case["constraints"], *case["golden"].values()]}
    result = driver.probe_case({"cpu": 2}, case, tmp_path / "probe")
    assert result["status"] == "analyzable_contexts_found"
    assert result["nomination_counts"]["selected_contexts"] == result["tree_contexts"] == result["accepted_contexts"] == 1
    assert result["effective_partition"]["effective_families"] == 3
    assert result["contexts"][0]["profile_accepted"]
    assert not result["accuracy_evaluated"] and not result["final_native_inference"]
    assert all(record(path) == pin for path, pin in before.items())
    proposal = json.loads(Path(result["contexts"][0]["proposal"]["path"]).read_text())
    assert proposal["joint_admission"]["baseline_high_pairs_lost"] == []
    assert proposal["joint_admission"]["novel_cross_high_pairs"]


def test_failed_context_is_retained_and_not_reclassified_as_inactivity(tmp_path, monkeypatch):
    case = case_fixture(tmp_path / "case")
    def failed(*args, **kwargs):
        raise RuntimeError("controlled local failure")
    monkeypatch.setattr(driver, "evaluate_local_context", failed)
    folder = tmp_path / "probe"
    with pytest.raises(RuntimeError, match="controlled local failure"):
        driver.probe_case({"cpu": 2}, case, folder)
    result = json.loads((folder / "probe.json").read_text())
    assert result["status"] == "failed" and result["contexts"][0]["status"] == "failed"
    assert result["error_type"] == "RuntimeError" and result["parent_resources"]["wall_s"] > 0
    with pytest.raises(FileExistsError):
        driver.probe_case({"cpu": 2}, case, folder)
