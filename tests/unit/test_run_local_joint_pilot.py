import json
import os
from pathlib import Path
import shutil

import numpy as np
import pytest

from benchmark_tools import run_local_joint_pilot as pilot
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from orthohmm.accuracy import write_accuracy_checkpoint
from orthohmm.local_rescue import nominate_local_contexts, evaluate_local_context, write_effective_candidates
from tests.unit.test_local_rescue_readback import accepted


@pytest.fixture
def prepared(accepted, tmp_path, monkeypatch):
    context, names, database, owners, sequences, tree_pin, files, config, constraints, pairs = accepted["arguments"]
    baseline = tmp_path / "baseline"
    shutil.copytree(accepted["baseline"].parent, baseline)
    for name in ("checkpoints", "gene_trees"):
        (baseline / "orthohmm_phylogeny" / name).mkdir(exist_ok=True)
    tree_pin = record(baseline / "orthohmm_phylogeny/species_tree.rooted.nwk")
    groups = [[g for g in names if g.endswith(f"family{i}")] for i in range(3)]
    groups.append([g for g in names if g.endswith("duplicate")])
    indices = {gene: i for i, gene in enumerate(names)}
    classes = {taxon: i for i, taxon in enumerate(sorted(set(owners.values())))}
    species = np.array([classes[owners[g]] for g in names], dtype=np.int32)
    hits = [(indices[a], indices[b], 2.) for left, right in ((groups[2], groups[3]), (groups[3], groups[2]))
            for a in left for b in right if owners[a] != owners[b]]
    q, t, scores = (np.array([row[i] for row in hits], dtype=np.float64 if i == 2 else np.int32) for i in range(3))
    checkpoint = write_accuracy_checkpoint(str(tmp_path / "initial"), names, species, q, t, scores)
    candidate = tmp_path / "original_candidates.txt"
    candidate.write_text("".join(" ".join(g) + "\n" for g in groups))
    constraint_path = tmp_path / "constraints.json"
    constraint_path.write_text("[]\n")
    nomination = nominate_local_contexts([[indices[g] for g in group] for group in groups], names, species, q, t, scores)
    assert len(nomination["contexts"]) == 1
    root = tmp_path / "probe/OrthoBench"
    root.mkdir(parents=True)
    nomination_path = root / "nomination.json"
    nomination_path.write_text(json.dumps(nomination))
    context, = nomination["contexts"]
    directory = root / "contexts" / context["context_id"]
    proposal = evaluate_local_context(context, names, database, owners, sequences, tree_pin["path"], files,
        config, directory, [], pairs, cpu=1)
    assert proposal["accepted"]
    effective = root / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"
    audit = write_effective_candidates(candidate, names, [proposal], effective, [])
    report = {"status": "analyzable_contexts_found", "dataset": "OrthoBench", "nomination": record(nomination_path),
        "contexts": [{"context_id": context["context_id"], "proposal": record(directory / "proposal.json"),
            "accepted": True, "tree_inferred": True, "profile_accepted": True}], "accepted_contexts": 1,
        "tree_contexts": 1, "candidate": record(effective), "effective_partition": audit}
    report_path = root / "probe.json"
    report_path.write_text(json.dumps(report))
    source = Path(pilot.__file__).resolve().parents[1]
    # This synthetic proposal was built above with current modules, not the historical probe.
    original_source = source
    case = {"dataset": "OrthoBench", "fasta_directory": str(accepted["root"] / "input"),
        "baseline_output": str(baseline), "checkpoint": str(checkpoint),
        "checkpoint_manifest": record(checkpoint / "manifest.json"),
        "fastas": [record(accepted["root"] / "input" / f) for f in files], "candidate": record(candidate),
        "constraints": record(constraint_path), "aligner": config.aligner, "tree_builder": config.tree_builder,
        "golden": {name: record(baseline / "orthohmm_phylogeny" / name) for name in (*pilot.PRODUCTS, "species_tree.rooted.nwk")},
        "expected": {"f1": .7, "precision": .8, "recall": .6}}
    output = tmp_path / "pilot"
    output.mkdir()
    pins = [record(p) for p in sorted(root.rglob("*")) if p.is_file()]
    plan = {"status": "prepared_local_joint_readback_and_native_unrun", "source_root": str(source), "repo": str(source),
        "output": str(output), "cpu": 16, "source_revision": "synthetic_fixture_not_frozen_biological_run", "case": case,
        "original_source": str(original_source), "case_report": record(report_path), "inputs": pins, "sources": [],
        "launch_ledger": pilot.committed(source, pilot.LEDGER)[1],
        "runtime": {"files": []}, "tool_path": os.environ["PATH"], "scope": "Actual16-gene integration, mock benchmark scores, not biological accuracy"}
    path = output / "execution_plan.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(pilot, "UNIVERSE", (16, 4))
    monkeypatch.setattr(pilot, "runtime", lambda: {"files": []})
    monkeypatch.setattr(pilot, "score_orthobench", lambda *a: {"scores": {"OrthoBench F1": .7}})
    monkeypatch.setattr(pilot, "readback_scores", lambda *a: ({"OrthoBench F1": .7}, {"precision": .8, "recall": .6}, []))
    return path, plan, report


def test_actual_nomination_profiles_cached_admission_and_final_inferred_native_pipeline(prepared):
    path, plan, _ = prepared
    result = pilot.run(path)
    assert result["status"] == "finished_local_readback_and_native_pilot_verified"
    off, actual = result["arms"]
    assert all(off["golden_matches"].values())
    assert off["status"] == actual["status"] == "independent_content_readback_verified"
    assert actual["accepted_contexts"] == ["LocalContext000001"]
    assert actual["native_pair_turnover"]["added_pairs"] > 0
    assert actual["final_effects"]["genes_with_changed_root_og_membership"] > 0
    assert not actual["advancement"]["eligible"]  # Mock neutral scores do not claim accuracy.
    assert actual["coverage"]["all_candidate_genes_preserved"]
    with pytest.raises(ValueError, match="fresh attempt"):
        pilot.run(path)


def test_failed_control_score_gate_prevents_local_native_scoring(prepared):
    path, plan, _ = prepared
    plan["case"]["expected"]["f1"] += .001
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="off control scores"):
        pilot.run(path)
    report = json.loads((Path(plan["output"]) / "results.json").read_text())
    assert report["status"] == "failed" and [row["label"] for row in report["arms"]] == ["off"]
    assert not (Path(plan["output"]) / "local_joint_v1").exists()


@pytest.mark.parametrize("kind", ["nomination", "context_list", "accepted", "partition"])
def test_complete_readback_rejects_modified_panel_or_partition(prepared, tmp_path, kind):
    _, plan, report = prepared
    if kind == "nomination":
        p = Path(report["nomination"]["path"])
        data = json.loads(p.read_text())
        data["counts"]["selected_contexts"] += 1
        p.write_text(json.dumps(data))
    elif kind == "context_list": report["contexts"] = []
    elif kind == "accepted": report["accepted_contexts"] = 2
    else: report["candidate"]["sha256"] = "wrong"
    Path(plan["case_report"]["path"]).write_text(json.dumps(report))
    with pytest.raises(ValueError):
        pilot.read_proposals(plan, tmp_path / "readback")


def test_final_effects_distinguishes_provisional_pairs_from_final_confidence(tmp_path):
    old = tmp_path / "baseline/orthohmm_working_res/orthohmm_edges_clustered.txt"
    old.parent.mkdir(parents=True)
    old.write_text("a b\nc d\n")
    new, confidence = tmp_path / "new.txt", tmp_path / "confidence.tsv"
    new.write_text("a c\nb\nd\n")
    confidence.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\tconfidence\na\tA\tc\tC\tmedium\n")
    contexts = [{"accepted": True, "context_id": "context", "genes": list("abcd"),
        "baseline_high_pairs": [["a", "b"], ["c", "d"]], "joint_admission": {"novel_cross_high_pairs": [["a", "c"]]}}]
    actual = pilot.final_effects({"baseline_output": str(old.parents[1])}, new, confidence, contexts)
    assert actual["genes_with_changed_root_og_membership"] == 4
    row, = actual["contexts"]
    assert row["baseline_high_pairs_absent_final"] == 2
    assert row["provisional_novel_pairs_present_final"] == 1 and row["provisional_novel_pairs_high_final"] == 0


def test_prepare_requires_frozen_reader_and_fresh_namespace(tmp_path):
    with pytest.raises(ValueError, match="frozen scientific reader"):
        pilot.prepare(tmp_path, tmp_path, tmp_path / "new")
    with pytest.raises(ValueError, match="fresh output"):
        pilot.prepare(Path("/repo"), Path("/source"), tmp_path)


def test_active_ledger_updates_do_not_invalidate_committed_launch(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *a, **kw: b'{"original": true}\n')
    path = tmp_path / pilot.LEDGER
    path.parent.mkdir(parents=True)
    path.write_text('{"current_status": "changed"}\n')
    value, binding = pilot.committed(tmp_path, pilot.LEDGER)
    assert value == {"original": True} and binding["revision"] == pilot.LAUNCH
    assert "git_blob_sha256" in binding and "bytes" not in binding


@pytest.mark.parametrize("state", ["COMPLETED", "RUNNING", "FAILED"])
def test_original_scheduler_state_is_checked_not_inferred_from_files(monkeypatch, state):
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *a, **kw:
        "JobIDRaw|State|ExitCode|AllocCPUS|ReqMem|NodeList\n24489|" + state + "|0:0|16|128G|bizon\n")
    if state == "COMPLETED": assert pilot.terminal("24489")["State"] == state
    else:
        with pytest.raises(ValueError, match="terminal successful"): pilot.terminal("24489")
