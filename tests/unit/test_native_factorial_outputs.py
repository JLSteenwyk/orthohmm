"""Native output semantics use existing tiny diagnostics, never rerun inference."""

from copy import deepcopy
import itertools
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from benchmark_tools import validate_native_factorial_outputs as review
from benchmark_tools.prepare_ob_candidate_neighborhood import record


FIXTURE = review.ROOT / "benchmark_tools/results/native_factorial_adapter_diagnostic_20261004"
EMBEDDED = str(Path(json.loads((FIXTURE / "p0_c0_r1/metrics.json").read_text())
                    ["metadata"]["fasta_directory"]).parent)
CELLS = ["p%d_c%d_r%d" % bits for bits in itertools.product(range(2), repeat=3)]


def context(cell, root=FIXTURE):
    metrics = json.loads((root / cell / "metrics.json").read_text())
    meta = metrics["metadata"]
    tools = json.loads((FIXTURE / "p0_c0_r1/metrics.json").read_text())["metadata"]
    return dict(cell=cell, output_root=str(root / cell), input_directory=meta["fasta_directory"],
        inputs=[record(p) for p in sorted(Path(meta["fasta_directory"]).iterdir())],
        native_order=["s0.fa", "s2.fa", "s1.fa", "s3.fa"], genes=26, proteomes=4,
        cpu=2, threads_per_worker=2, command=metrics["command"], cwd=metrics["cwd"],
        aligner=tools["aligner"], tree_builder=tools["tree_builder"])


@pytest.mark.parametrize("cell", CELLS)
def test_retained_actual_native_cells(cell, tmp_path):
    result = review.validate_semantics(copy_context(tmp_path, cell))
    assert result["native_outputs_validated"] is True
    assert result["input_genes"] == 26
    assert result["orthogroups"] == 6
    assert result["checkpoint"]["hits"] == 116
    assert result["checkpoint"]["order"] == "lexical_gene_ids"
    assert result["accuracy_evaluated"] is False
    assert result["resource_measurements_admitted"] is False
    assert result["next_identity_authorized"] is False
    if cell.endswith("r1"):
        assert result["phylogeny"]["native_pair_rows"] == 43
    else:
        assert result["phylogeny"] is None
    for pin in result["checked_files"]:
        assert record(pin["path"]) == pin


def copy_context(tmp_path, cell):
    shutil.copytree(FIXTURE / "fixture", tmp_path / "fixture")
    shutil.copytree(FIXTURE / cell, tmp_path / cell)
    # Relocate only explicit embedded artifact paths in this synthetic copy.
    def relocate(value):
        if isinstance(value, str) and value.startswith(EMBEDDED + "/"):
            return str(tmp_path) + value[len(EMBEDDED):]
        if isinstance(value, list):
            return [relocate(v) for v in value]
        if isinstance(value, dict):
            return {k: relocate(v) for k, v in value.items()}
        return value
    for path in (tmp_path / cell).rglob("*.json"):
        value = relocate(json.loads(path.read_text()))
        path.write_text(json.dumps(value))
    ctx = context(cell, tmp_path)
    return ctx


@pytest.fixture
def copied(tmp_path):
    ctx = copy_context(tmp_path, "p1_c1_r1")
    review.validate_semantics(ctx)
    return ctx


def change_json(path, mutate):
    value = json.loads(path.read_text())
    mutate(value)
    path.write_text(json.dumps(value))


@pytest.mark.parametrize("key,value", [
    ("cpu_budget", 1), ("search_kmer_k", 3), ("search_max_candidates_per_query", 50),
    ("evalue_threshold", .001), ("cpm_resolution", .2), ("leiden_seed", 5),
    ("phylogeny", "off"), ("species_tree_mode", "supplied"), ("species_tree_rooting", "midpoint"),
    ("phylogeny_pair_rule", "strict"), ("fasta_directory", "/tmp/other"),
])
def test_frozen_metadata_tampering(copied, key, value):
    path = Path(copied["output_root"]) / "metrics.json"
    change_json(path, lambda m: m["metadata"].update({key: value}))
    with pytest.raises(ValueError, match="frozen settings"):
        review.validate_semantics(copied)


@pytest.mark.parametrize("field", ["command", "cwd", "status"])
def test_completion_identity_tampering(copied, field):
    path = Path(copied["output_root"]) / "metrics.json"
    change_json(path, lambda m: m.update({field: "wrong"}))
    with pytest.raises(ValueError, match="completion"):
        review.validate_semantics(copied)


def test_missing_stage(copied):
    change_json(Path(copied["output_root"]) / "metrics.json", lambda m: m["stages"].pop("search"))
    with pytest.raises(ValueError, match="stage set"):
        review.validate_semantics(copied)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1, 0, True])
def test_invalid_stage_duration(copied, value):
    change_json(Path(copied["output_root"]) / "metrics.json",
                lambda m: m["stages"]["search"].update(wall_s=value))
    with pytest.raises(ValueError, match="stage duration"):
        review.validate_semantics(copied)


def test_native_factor_tampering(copied):
    change_json(Path(copied["output_root"]) / "metrics.json",
                lambda m: m["metadata"]["native_factorial"].update(production_globals_modified=True))
    with pytest.raises(ValueError, match="frozen settings"):
        review.validate_semantics(copied)


@pytest.mark.parametrize("name", ["genes", "species", "orthogroups", "phylogeny_ortholog_pairs"])
def test_count_tampering(copied, name):
    change_json(Path(copied["output_root"]) / "metrics.json",
                lambda m: m["counts"].update({name: m["counts"][name] + 1}))
    with pytest.raises(ValueError):
        review.validate_semantics(copied)


def test_boolean_count_rejected(copied):
    change_json(Path(copied["output_root"]) / "metrics.json", lambda m: m["counts"].update(genes=True))
    with pytest.raises(ValueError, match="integer count"):
        review.validate_semantics(copied)


@pytest.mark.parametrize("mutation", ["duplicate", "missing", "unknown"])
def test_materialized_universe_tampering(copied, mutation):
    path = Path(copied["output_root"]) / "native/orthohmm_orthogroups.txt"
    lines = path.read_text().splitlines()
    gene = lines[0].split()[1]
    if mutation == "duplicate":
        lines[1] += " " + gene
    elif mutation == "missing":
        lines[0] = lines[0].replace(" " + gene, "", 1)
    else:
        lines[0] += " not_an_input_gene"
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        review.validate_semantics(copied)


def reseal_checkpoint(path, name):
    pin = record(path / name)
    change_json(path / "manifest.json",
                lambda m: m["files"].update({name: {k: pin[k] for k in ("bytes", "sha256")}}))


@pytest.mark.parametrize("name,mutation", [
    ("hit_scores.npy", "nan"), ("hit_queries.npy", "negative"),
    ("hit_targets.npy", "out_of_range"), ("gene_to_species.npy", "ownership"),
    ("hit_queries.npy", "dtype"), ("hit_scores.npy", "shape"),
])
def test_semantically_bad_resealed_arrays(copied, name, mutation):
    path = Path(copied["output_root"]) / "native/orthohmm_working_res/high_sensitivity_checkpoint"
    array = np.load(path / name, allow_pickle=False)
    if mutation == "nan":
        array[0] = float("nan")
    elif mutation == "negative":
        array[0] = -1
    elif mutation == "out_of_range":
        array[0] = copied["genes"]
    elif mutation == "ownership":
        array[0] = 100
    elif mutation == "dtype":
        array = array.astype("int64")
    else:
        array = array.reshape(1, -1)
    np.save(path / name, array, allow_pickle=False)
    reseal_checkpoint(path, name)
    with pytest.raises(ValueError, match="Checkpoint|checkpoint"):
        review.validate_semantics(copied)


def test_resealed_wrong_checkpoint_order(copied):
    path = Path(copied["output_root"]) / "native/orthohmm_working_res/high_sensitivity_checkpoint"
    names = (path / "gene_names.txt").read_text().splitlines()
    (path / "gene_names.txt").write_text("\n".join(reversed(names)) + "\n")
    reseal_checkpoint(path, "gene_names.txt")
    with pytest.raises(ValueError, match="lexical"):
        review.validate_semantics(copied)


def test_checkpoint_file_inventory(copied):
    path = Path(copied["output_root"]) / "native/orthohmm_working_res/high_sensitivity_checkpoint/manifest.json"
    change_json(path, lambda m: m["files"].pop("hit_scores.npy"))
    with pytest.raises(ValueError, match="inventory"):
        review.validate_semantics(copied)


def test_checkpoint_bounded_chunks(copied):
    evidence = review.Evidence()
    owners, _, _ = review.universe(copied, evidence)
    root = Path(copied["output_root"])
    metrics = json.loads((root / "metrics.json").read_text())
    assert review.checkpoint(root / "native/orthohmm_working_res/high_sensitivity_checkpoint",
                             owners, metrics, evidence, chunk_size=3)["hits"] == 116
    evidence.finish()


@pytest.mark.parametrize("mutation", ["duplicate", "reverse", "species", "unknown"])
def test_native_pair_tampering(copied, mutation):
    path = Path(copied["output_root"]) / "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
    lines = path.read_text().splitlines()
    row = lines[1].split("\t")
    if mutation == "duplicate":
        lines.insert(2, lines[1])
    elif mutation == "reverse":
        lines[1] = "\t".join(row[2:] + row[:2])
    elif mutation == "species":
        row[1] = row[3]
        lines[1] = "\t".join(row)
    else:
        row[0] = "absent_gene"
        lines[1] = "\t".join(row)
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        review.validate_semantics(copied)


def test_resealed_wrong_tree_taxa(copied):
    directory = Path(copied["output_root"]) / "native/orthohmm_phylogeny"
    tree = directory / "species_tree.rooted.nwk"
    tree.write_text("(s0:1,s1:1,s2:1,wrong:1);\n")
    change_json(directory / "provenance_manifest.json",
                lambda m: m.update(species_tree_sha256=record(tree)["sha256"]))
    with pytest.raises(ValueError, match="taxon coverage"):
        review.validate_semantics(copied)


def test_candidate_parameter_tampering(copied):
    change_json(Path(copied["output_root"]) / "metrics.json",
                lambda m: m["metadata"]["phylogeny_candidate_profile"]["parameters"].update(max_iterations=3))
    with pytest.raises(ValueError, match="Candidate settings"):
        review.validate_semantics(copied)


def test_seed_sidecar_duplicate(copied):
    path = Path(copied["output_root"]) / "native/orthohmm_working_res/phylogeny_candidate_seeds.tsv"
    path.write_text(path.read_text().replace("Seed0000001", "Seed0000000"))
    with pytest.raises(ValueError, match="seed coverage"):
        review.validate_semantics(copied)


def test_phylogeny_checkpoint_reuse(copied):
    path = Path(copied["output_root"]) / "metrics.json"
    change_json(path, lambda m: m["counts"].update(phylogeny_species_tree_checkpoint_hit=True))
    with pytest.raises(ValueError, match="checkpoint reused"):
        review.validate_semantics(copied)


def test_symlink_evidence_rejected(copied, tmp_path):
    path = Path(copied["output_root"]) / "native/orthohmm_orthogroups.txt"
    moved = tmp_path / "elsewhere.txt"
    path.rename(moved)
    path.symlink_to(moved)
    with pytest.raises(ValueError, match="direct absolute"):
        review.validate_semantics(copied)


def test_evidence_mutation_after_binding(tmp_path):
    path = tmp_path / "evidence.json"
    path.write_text("{}")
    evidence = review.Evidence()
    evidence.bind(path)
    path.write_text("{\"changed\":true}")
    with pytest.raises(ValueError):
        evidence.finish()


@pytest.mark.parametrize("state,exit_code", [("RUNNING", "0:0"), ("FAILED", "1:0"), ("COMPLETED", "1:0")])
def test_production_requires_clean_terminal(monkeypatch, state, exit_code):
    request = dict(plan={"path": "/plan"}, index=0, job_id=22427)
    monkeypatch.setattr(review, "read", lambda ref: request if ref == {"path": "/request"} else {})
    monkeypatch.setattr(review, "validate_plan", lambda plan: [])
    monkeypatch.setattr(review, "validate_request", lambda *args: None)
    monkeypatch.setattr(review, "verify_terminal", lambda job: dict(verified=dict(fields=dict(JobState=state, ExitCode=exit_code))))
    with pytest.raises(ValueError, match="successful terminal"):
        review.validate({"path": "/request"})


def test_input_duplicate_ownership(copied):
    ctx = deepcopy(copied)
    path = Path(ctx["input_directory"]) / "s0.fa"
    path.write_text(path.read_text() + ">s0_f0\nACDE\n")
    ctx["inputs"] = [record(Path(ctx["input_directory"]) / Path(ref["path"]).name) for ref in ctx["inputs"]]
    with pytest.raises(ValueError, match="duplicate input gene"):
        review.validate_semantics(ctx)


def test_nonzero_candidate_merge_kernel(tmp_path):
    working = tmp_path / "working"
    working.mkdir()
    partitions = working / "phylogeny_candidate_superfamilies.txt"
    partitions.write_text("a b c\nd\n")
    seeds = working / "phylogeny_candidate_seeds.tsv"
    seeds.write_text("candidate_family\tseed_families\nFamily0000000\tSeed0000000,Seed0000001\nFamily0000001\tSeed0000002\n")
    trace = working / "phylogeny_candidate_merges.json"
    trace.write_text(json.dumps([dict(source_genes=["a"], target_genes=["b", "c"])]))
    profile = dict(profile="satellite_v2", parameters=review.SATELLITE,
        membership_policy="high_confidence_pair", candidate_families=2, seed_families=3,
        merges=1, iterations=1, candidate_checkpoint=str(partitions), seed_sidecar=str(seeds),
        merge_trace_sidecar=str(trace))
    counts = dict(phylogeny_seed_families=3, phylogeny_candidate_merges=1)
    owners = dict(a="s1", b="s2", c="s3", d="s1")
    evidence = review.Evidence()
    path, groups, constraints = review.candidate(dict(phylogeny_candidate_profile=profile),
        counts, working, owners, evidence)
    assert path == partitions and len(groups) == 2 and constraints == 1
    evidence.finish()
    trace.write_text(json.dumps([dict(source_genes=["a"], target_genes=["d"])]))
    with pytest.raises(ValueError, match="crosses candidate"):
        review.candidate(dict(phylogeny_candidate_profile=profile), counts, working, owners, review.Evidence())


@pytest.mark.parametrize("cell", ["p0_c0_r0", "p0_c1_r0"])
def test_r_off_cannot_hide_reconciliation_metadata(cell, tmp_path):
    ctx = copy_context(tmp_path, cell)
    change_json(Path(ctx["output_root"]) / "metrics.json", lambda m: m["metadata"].update(phylogeny="reconcile"))
    with pytest.raises(ValueError, match="R-off"):
        review.validate_semantics(ctx)


@pytest.mark.parametrize("representation", ["controller", "accounting"])
def test_production_context_and_receipt_binding(copied, tmp_path, monkeypatch, representation):
    semantic = review.validate_semantics(copied)
    run = dict(copied, index=0, input_directory=str(tmp_path / "source"))
    baseline = dict(core_root="/frozen/core", tool_entrypoints={
        "orthohmm_python": dict(absolute_path="/private/python"),
        "mafft": dict(absolute_path=copied["aligner"]),
        "FastTree": dict(absolute_path=copied["tree_builder"])})
    baseline_file = tmp_path / "baseline.json"
    baseline_file.write_text(json.dumps(baseline))
    plan_file = tmp_path / "plan.json"
    plan_file.write_text(json.dumps(dict(baseline=record(baseline_file), helper_sources=[])))
    plan_ref = record(plan_file)
    request_file = tmp_path / "request.json"
    request_file.write_text(json.dumps(dict(plan=plan_ref, index=0, job_id=22427)))
    preparation = dict(status="fresh_factorial_inputs_prepared", genes=26,
        gene_ownership_sha256=semantic["gene_ownership_sha256"], per_species_counts=semantic["per_species_counts"])
    execution = dict(status="native_factorial_completed_pending_output_review", plan=plan_ref,
        index=0, cell=copied["cell"], factors=semantic["factors"],
        native_order=copied["native_order"], automatic_retry=False)
    root = Path(copied["output_root"])
    (root / "preparation.json").write_text(json.dumps(preparation))
    (root / "native_execution.json").write_text(json.dumps(execution))
    monkeypatch.setattr(review, "validate_plan", lambda plan: [run])
    monkeypatch.setattr(review, "validate_request", lambda *args: None)
    verified = dict(fields=dict(JobState="COMPLETED", ExitCode="0:0")) if representation == "controller" else dict(State="COMPLETED", ExitCode="0:0")
    monkeypatch.setattr(review, "verify_terminal", lambda job: dict(verified=verified))
    observed = []
    monkeypatch.setattr(review, "validate_semantics", lambda ctx: observed.append(ctx) or semantic)
    result = review.validate(record(request_file))
    assert result["terminal_scheduler_confirmed"] is True
    assert result["terminal_reviewed"] is False
    assert result["next_identity_authorized"] is False
    assert observed[0]["cpu"] == 32 and observed[0]["threads_per_worker"] == 4
    assert observed[0]["command"] == ["/private/python", str(review.ROOT / "benchmark_tools/run_native_factorial_cost.py"),
        "--native", "--plan", plan_ref["path"], "--plan-sha256", plan_ref["sha256"], "--index", "0"]
    assert observed[0]["cwd"] == "/frozen/core"
    change_json(root / "native_execution.json", lambda e: e.update(automatic_retry=True))
    with pytest.raises(ValueError, match="execution receipt"):
        review.validate(record(request_file))
