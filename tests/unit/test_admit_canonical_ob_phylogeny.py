import json

import pytest

from benchmark_tools.admit_canonical_ob_phylogeny import verify_artifacts
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_canonical_ob_phylogeny import replay_arguments
from benchmark_tools.verify_ygob_validation import require_completed_job


def write(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def evidence(tmp_path):
    inputs = tmp_path / "input"
    inputs.mkdir()
    (inputs / "one.fa").write_text(">a\nACDE\n")
    candidate, constraints = tmp_path / "candidate", tmp_path / "constraints"
    candidate.write_text("a\n")
    constraints.write_text("[]")
    source = tmp_path / "driver.py"
    source.write_text("# fixture\n")
    plan = dict(output=str(tmp_path), inference=str(tmp_path / "inference"), inputs=str(inputs),
        attempts=1, checkpoint_reuse=False, scoring=False, cpu=1, python="/python",
        aligner="/mafft", tree_builder="/FastTree", runtime={"test": "runtime"},
        source=record(source), candidate=record(candidate), constraints=record(constraints),
        checked_records=[record(p) for p in (inputs / "one.fa", candidate, constraints)])
    write(tmp_path / "plan.json", plan)
    pinned = record(tmp_path / "plan.json")
    command = [str(source), "--run", pinned["path"], "--plan-sha256", pinned["sha256"]]
    write(tmp_path / "started.json", dict(plan=pinned, runtime=plan["runtime"],
        command=command, arguments=replay_arguments(plan)))
    phylo = tmp_path / "inference/orthohmm_phylogeny"
    phylo.mkdir(parents=True)
    summary = dict(checkpoint_hits=0, remapped_checkpoint_hits=0, species_tree_checkpoint_hit=False,
                   species_tree_families=1, reconciled_families=1)
    rules = dict(species_tree_mode="infer", species_tree_rooting="min_variance",
                 root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    manifest = dict(input_cluster_sha256=plan["candidate"]["sha256"], results=summary,
                    species_tree_taxa=["one"], membership_reconciliation={"policy": "high_confidence_pair"}, **rules)
    write(phylo / "provenance_manifest.json", manifest)
    write(phylo / "reconciliation_summary.json", summary)
    (phylo / "orthohmm_root_hogs.tsv").write_text("fixture only; gate does not validate groups\n")
    replay = dict(status="complete", command=[plan["python"], *command],
        input=dict(candidate_clusters=plan["candidate"], membership_constraints=plan["constraints"],
                   fasta_directory=str(inputs), files=["one.fa"]),
        outputs={k:record(phylo / n) for k,n in (("manifest","provenance_manifest.json"),
            ("summary","reconciliation_summary.json"),("root_hogs","orthohmm_root_hogs.tsv"))},
        counts=summary, parameters=dict(checkpoint_source=None, species_tree=None,
            explicit_unconstrained_ablation=False, **rules))
    write(tmp_path / "replay.json", replay)
    write(tmp_path / "native_complete.json", dict(plan=pinned, accuracy_evaluated=False,
        status="fresh_canonical_phylogeny_complete_pending_readback", replay=record(tmp_path / "replay.json")))
    return tmp_path, pinned["sha256"]


def test_binding_is_not_scientific_admission(evidence):
    result = verify_artifacts(*evidence, 1, 1)
    assert result["genes"] == 1 and result["scientific_scores_admitted"] is False


@pytest.mark.parametrize("change", ["plan", "failure", "runtime", "arguments", "command",
                                    "completion", "replay", "output", "input", "extra_input"])
def test_mutation_rejected(evidence, change):
    directory, sha = evidence
    if change == "failure":
        write(directory / "failure.json", {})
    elif change in ("input", "extra_input"):
        name = "one.fa" if change == "input" else "extra.fa"
        (directory / "input" / name).write_text(">a\nAAAA\n")
    elif change == "output":
        (directory / "inference/orthohmm_phylogeny/orthohmm_root_hogs.tsv").write_text("changed")
    else:
        file, key = {"plan": ("plan.json", "cpu"), "runtime": ("started.json", "runtime"),
            "arguments": ("started.json", "arguments"), "command": ("started.json", "command"),
            "completion": ("native_complete.json", "accuracy_evaluated"),
            "replay": ("replay.json", "status")}[change]
        path = directory / file
        value = json.loads(path.read_text())
        value[key] = "changed"
        write(path, value)
    with pytest.raises(ValueError):
        verify_artifacts(directory, sha, 1, 1)


@pytest.mark.parametrize("change", ["checkpoint", "remapped", "species_checkpoint", "supplied_tree",
                                    "root_rule", "pair_rule", "constraint_policy", "candidate"])
def test_rehashed_semantic_mutation_rejected(evidence, change):
    directory, sha = evidence
    phylo = directory / "inference/orthohmm_phylogeny"
    replay = json.loads((directory / "replay.json").read_text())
    manifest = json.loads((phylo / "provenance_manifest.json").read_text())
    summary = manifest["results"]
    if change in ("checkpoint", "remapped", "species_checkpoint"):
        key = {"checkpoint": "checkpoint_hits", "remapped": "remapped_checkpoint_hits",
               "species_checkpoint": "species_tree_checkpoint_hit"}[change]
        summary[key] = True if change == "species_checkpoint" else 1
        replay["counts"] = summary
    elif change == "supplied_tree":
        replay["parameters"]["species_tree"] = "/external.nwk"
    elif change in ("root_rule", "pair_rule"):
        key = "root_duplication_rule" if change == "root_rule" else "pair_orthology_rule"
        manifest[key] = replay["parameters"][key] = "changed"
    elif change == "constraint_policy":
        manifest["membership_reconciliation"]["policy"] = "none"
    else:
        manifest["input_cluster_sha256"] = "0" * 64
    write(phylo / "provenance_manifest.json", manifest)
    write(phylo / "reconciliation_summary.json", summary)
    for key, name in (("manifest", "provenance_manifest.json"), ("summary", "reconciliation_summary.json")):
        replay["outputs"][key] = record(phylo / name)
    write(directory / "replay.json", replay)
    completion = json.loads((directory / "native_complete.json").read_text())
    completion["replay"] = record(directory / "replay.json")
    write(directory / "native_complete.json", completion)
    with pytest.raises(ValueError):
        verify_artifacts(directory, sha, 1, 1)


@pytest.mark.parametrize("state,code", [("RUNNING", "0:0"), ("FAILED", "1:0"), ("COMPLETED", "1:0")])
def test_scheduler_must_be_successful_terminal(state, code):
    with pytest.raises(ValueError):
        require_completed_job(f"JobIDRaw|State|ExitCode|Elapsed\n22324|{state}|{code}|00:01:00\n", 22324)
