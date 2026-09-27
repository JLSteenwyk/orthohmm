import copy

import pytest

from benchmark_tools.readback_qfo_canonical_phylogeny import (
    config_arguments, record, verify_cache, verify_execution,
)
from benchmark_tools.qfo_phylogeny_cache import cache_paths


def cache_fixture(tmp_path):
    source, target = tmp_path / "source", tmp_path / "target"
    family = "Family0000001"
    files = []
    for relative in cache_paths(family):
        path = source / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("bytes")
        files.append(dict(relative=relative, source=record(path)))
    tool = tmp_path / "mafft/bin/mafft"
    tool.parent.mkdir(parents=True)
    tool.write_text("tool")
    tool_record = record(tool)
    plan = dict(aligner=str(tool), tree_builder=str(tool), tool_records=[tool_record], checked_records=[tool_record])
    gate = dict(source_directory=str(source), admitted_outputs=[r["source"] for r in files])
    cache = dict(source_directory=str(source), tool_records=[tool_record], files=files, included_families=[family])
    copies = [dict(r["source"], path=str(target / r["relative"])) for r in files]
    return plan, gate, cache, copies, target


def test_cache_copy_pre_execution_receipt(tmp_path):
    # Final checkpoint bytes may differ; this verifies their pre-execution copies.
    verify_cache(*cache_fixture(tmp_path))


@pytest.mark.parametrize("mutation", ["tool", "source", "duplicate_family", "extra_file", "copy", "unadmitted"])
def test_reject_cache_receipt_changes(tmp_path, mutation):
    args = list(cache_fixture(tmp_path))
    if mutation == "tool":
        args[2]["tool_records"] = []
    elif mutation == "source":
        args[2]["source_directory"] = "/other"
    elif mutation == "duplicate_family":
        args[2]["included_families"] *= 2
    elif mutation == "extra_file":
        args[2]["files"].append(args[2]["files"][0])
    elif mutation == "copy":
        args[3][0]["sha256"] = "changed"
    else:
        args[1]["admitted_outputs"] = []
    with pytest.raises(ValueError):
        verify_cache(*args)


def execution_fixture():
    source = dict(path="/repo/benchmark_tools/run_qfo_canonical_phylogeny.py", bytes=1, sha256="source")
    identity = dict(path="/plan", bytes=1, sha256="plan")
    receipt = dict(path="/started", bytes=1, sha256="start")
    cache_record, copies_record = dict(path="/cache"), dict(path="/copies")
    plan = dict(repo="/repo", checked_records=[source], python="/python", environment={},
                aligner="/mafft", tree_builder="/FastTree", reuse_admission=dict(path="/admission"))
    started = dict(plan=identity, job_id="99", source=source, executable=plan["python"],
        environment={}, config=config_arguments(plan), reuse_admission=plan["reuse_admission"], attempts=1)
    summary = dict(checkpoint_hits=1, remapped_checkpoint_hits=0,
                   species_tree_checkpoint_hit=False, species_tree_families=2)
    complete = dict(status="canonical_phylogeny_complete_pending_readback", started=receipt,
        plan=identity, cache=cache_record, copies=copies_record, accuracy_evaluated=False, summary=summary)
    cache = dict(included_families=["Family0000001"], species_tree_cache_copied=False,
                 reconciliation_outputs_copied=False)
    return plan, identity, 99, started, receipt, complete, cache, cache_record, copies_record


def test_bound_canonical_execution():
    verify_execution(*execution_fixture())


@pytest.mark.parametrize("index,key,value", [(3,"plan",{}), (3,"job_id","98"),
    (3,"attempts",2), (3,"environment",{"changed":"1"}), (3,"config",{}),
    (3,"reuse_admission",{}), (5,"started",{}), (5,"copies",{}), (5,"accuracy_evaluated",True),
    (6,"species_tree_cache_copied",True), (6,"reconciliation_outputs_copied",True)])
def test_canonical_binding_mutations(index, key, value):
    args = copy.deepcopy(execution_fixture())
    args[index][key] = value
    with pytest.raises(ValueError):
        verify_execution(*args)


@pytest.mark.parametrize("key,value", [("checkpoint_hits",2), ("remapped_checkpoint_hits",1),
    ("species_tree_checkpoint_hit",True), ("species_tree_families",0)])
def test_wrong_reuse_counts(key, value):
    args = execution_fixture()
    args[5]["summary"][key] = value
    with pytest.raises(ValueError):
        verify_execution(*args)
