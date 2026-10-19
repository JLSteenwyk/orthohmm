from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import statistics

import pytest

from benchmark_tools import report_projection_followup as reader
from benchmark_tools import run_projection_followup as followup
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.report_hillclimb import METRIC_NAMES, endpoint_vector
from benchmark_tools.run_membership_hillclimb import write_json
from tests.unit.test_profile_evidence_readback import arm_fixture, timing
from tests.unit.test_projection_followup_plan import evidence


def chosen(passing=None, inactive=()):
    report, producer = evidence(passing, inactive)
    return reader.planner.selection(report, producer)


def panel_fixture(tmp_path, monkeypatch, dataset="ThreeKingdoms", passing=None):
    selected = chosen(passing)
    rows, natives, expectations = [], {}, {}
    plan = None
    for label in selected["execution_arms"]:
        policy, mode = reader.recovery_configuration(label)
        activation = {"arm": label, "projection_policy": policy, "candidate_digest": "fixed_candidates",
            "completely_measured_queries": 3, "applicability": "completely_measured"}
        old, row, native, _, folder = arm_fixture(tmp_path, monkeypatch, label, experiment="projection_recovery", activation=activation)
        tree = folder / "native/orthohmm_phylogeny/species_tree.rooted.nwk"
        tree.write_text("(a,b);\n")
        seed = folder / "graph/orthogroups_profiles_refined.txt"
        seed.parent.mkdir()
        seed.write_text("a b\n")
        row["candidate"]["seed_partition"] = record(seed)
        graph_path = folder / "graph_metrics.json"
        graph = json.loads(graph_path.read_text())
        graph["stages"][0]["output"] = record(seed)
        write_json(graph_path, graph)
        score = .7 if label == "legacy_off" else .701
        scores = {metric: score for metric in METRIC_NAMES[dataset].values()}
        details = {"fixture": "mocked native inference/scoring; real partitions/commands/receipts"}
        if dataset == "QfO":
            details = {metric: {"precision": .8, "recall": .6} if metric.endswith(" F") else
                {"assessed_relations": 10, "similarity": score} for metric in scores}
        case = {"dataset": dataset, "case": old["case"], "checkpoint": old["case"]["checkpoint"],
            "checkpoint_sha": old["case"]["checkpoint_manifest"]["sha256"], "candidate": row["candidate"],
            "golden": {**old["case"]["golden"], "species_tree.rooted.nwk": record(tree)},
            "control": {"scores": endpoint_vector(dataset, scores)}}
        if plan is None:
            plan = {"dataset": dataset, "scientific_source": old["source_root"], "scientific_revision": followup.SCIENCE,
                "output": old["output"], "cpu": 16, "case": case, "selection": selected, "arms": selected["execution_arms"],
                "hmmer_tools": {}, "inputs": [], "scoring_ids": {name: followup.scoring_id(Path(old["output"]), name)
                    for name in selected["execution_arms"] if name not in ("legacy_off", "projection_hmmer_score_only")}}
        graph_command, native_command = reader.arm_commands(plan, folder, label)
        write_json(folder / "preflight.json", {"graph": graph_command, "native": native_command})
        row.update(dataset=dataset, role=selected["roles"][label], status=followup.SUCCESS,
            graph_execution=timing(folder, "graph", graph_command), native_execution=timing(folder, "native", native_command),
            graph_report=record(graph_path), graph_stage_signature=[["fixture", seed.stat().st_size, record(seed)["sha256"]]],
            scores=scores, details=details, metric_vector=endpoint_vector(dataset, scores), checked=[],
            rooted_species_tree=record(tree), native_pair_turnover={"added_pairs": 0, "files": [], "wall_s": 1.},
            golden_matches={"product.tsv": True, "species_tree.rooted.nwk": True})
        if dataset == "QfO": row.update(conversion={"fixture": True}, metric_files=[], references=[])
        write_json(folder / "assessment_receipt.json", {key: row[key] for key in ("scores", "details", "metric_vector", "checked")})
        row["checked"].append(record(folder / "assessment_receipt.json"))
        if label == "projection_hmmer_score_only":
            row["measurement_scoring_reuse"] = record(Path(plan["output"]) / "projection_off/assessment_receipt.json")
        native.update(genes=followup.UNIVERSES[dataset][0], species=followup.UNIVERSES[dataset][1],
            products={"product.tsv": record(folder / "native/orthohmm_phylogeny/product.tsv")})
        natives[label] = native
        expectations[label] = (deepcopy(scores), deepcopy(details), deepcopy(activation))
        rows.append(row)
    if dataset == "QfO":
        root = tmp_path / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
        root.mkdir(parents=True)
        write_json(root / "execution_plan.json", {"scoring_id": "canonical_fixture"})
        original = {"status": "canonical_qfo_scoring_recovered_independently_verified", "execution_plan": record(root / "execution_plan.json"),
            **{key: rows[0][key] for key in ("scores", "conversion", "metric_files", "references")}}
        write_json(root / "results.json", original)
        rows[0]["control_scoring_reuse"] = record(root / "results.json")
        plan["inputs"] += [record(root / name) for name in ("results.json", "execution_plan.json")]
    calls = []
    monkeypatch.setattr(reader, "replay_native", lambda row, *args: natives[row["label"]])
    monkeypatch.setattr(reader, "validate_activation", lambda graph, label, *args, **kwargs: expectations[label][2])
    monkeypatch.setattr(reader, "native_pair_turnover", lambda *args: {"added_pairs": 0, "files": [], "wall_s": 2.})
    def assess(repo, observed_dataset, row, native, folder, participant):
        assert observed_dataset == dataset
        calls.append((row["label"], folder, participant))
        return deepcopy(expectations[row["label"]][0]), deepcopy(expectations[row["label"]][1]), []
    monkeypatch.setattr(reader, "readback_scores", assess)
    monkeypatch.setattr(reader, "functional_populations", lambda folder, *args: ({"FAS": {"fixture_folder": str(folder)}}, []))
    return plan, rows, calls


def admitted_panel(tmp_path, monkeypatch, **kwargs):
    plan, rows, calls = panel_fixture(tmp_path, monkeypatch, **kwargs)
    cells = []
    for row in rows:
        projection = next((cell for cell in cells if cell["label"] == "projection_off"), None)
        cells.append(reader.admit_arm(tmp_path, plan, row, ["a", "b"], projection))
    return plan, rows, cells, calls


@pytest.mark.parametrize("dataset", reader.DATASETS)
def test_all_five_arms_are_independently_checked_with_exact_measurement_reuse(tmp_path, monkeypatch, dataset):
    plan, rows, cells, calls = admitted_panel(tmp_path, monkeypatch, dataset=dataset)
    deltas = reader.verify_dataset(cells, plan["selection"])
    assert len(deltas) == 4 and all(cell["checked"] for cell in cells)
    projection, measurement = cells[1], cells[3]
    assert measurement["scores"] == projection["scores"] and measurement["details"] == projection["details"]
    assert measurement["assessment_reuse"] == projection["assessment_receipt"]
    assert calls[3][1] == Path(plan["output"]) / "projection_off"
    assert cells[0]["scores"] == plan["case"]["control"]["scores"]
    if dataset == "QfO":
        assert calls[0][2] == "ohmm_hillclimb_r1_canonical_fixture"
        assert calls[3][2] == "ohmm_hillclimb_r1_" + plan["scoring_ids"]["projection_off"]
    else: assert all(call[2] is None for call in calls)


@pytest.mark.parametrize("mutation", ["command", "activation", "graph", "signature", "candidate", "seed", "coverage", "tree", "summary", "universe", "scores", "details", "vector", "turnover", "receipt", "role", "reuse"])
def test_mutated_arm_evidence_cannot_be_admitted(tmp_path, monkeypatch, mutation):
    plan, rows, _ = panel_fixture(tmp_path, monkeypatch)
    row = rows[2]
    folder = Path(plan["output"]) / row["label"]
    if mutation == "command": row["graph_execution"]["command"].append("--unexpected")
    elif mutation == "activation": row["activation"] = {"changed": True}
    elif mutation == "graph": row["graph_report"]["sha256"] = "wrong"
    elif mutation == "signature": row["graph_stage_signature"][0][2] = "wrong"
    elif mutation == "candidate": row["candidate"]["candidate_partition"]["sha256"] = "wrong"
    elif mutation == "seed": row["candidate"]["seed_partition"] = record(folder / "graph_metrics.json")
    elif mutation == "coverage": row["coverage"]["all_candidate_genes_preserved"] = False
    elif mutation == "tree": row["rooted_species_tree"]["sha256"] = "wrong"
    elif mutation == "summary": row["native_summary"] = {"wrong": True}
    elif mutation == "universe":
        original = reader.replay_native
        monkeypatch.setattr(reader, "replay_native", lambda *a: {**original(*a), "genes": 1})
    elif mutation == "scores": row["scores"]["Three Kingdoms F1"] += .01
    elif mutation == "details": row["details"]["unexpected"] = True
    elif mutation == "vector": row["metric_vector"][reader.ENDPOINTS[-1]] += .01
    elif mutation == "turnover": row["native_pair_turnover"]["added_pairs"] = 1
    elif mutation == "receipt": write_json(folder / "assessment_receipt.json", {"changed": True})
    elif mutation == "role": row["role"] = "advancing_candidate_changed"
    else: row["control_scoring_reuse"] = {"unexpected": True}
    with pytest.raises(ValueError): reader.admit_arm(tmp_path, plan, row, ["a", "b"])


@pytest.mark.parametrize("mutation", ["no_projection", "copied_points", "reuse_pin", "resampled", "native", "constraints"])
def test_measurement_requires_actual_repaired_assessment_and_native_identity(tmp_path, monkeypatch, mutation):
    plan, rows, _ = panel_fixture(tmp_path, monkeypatch, dataset="QfO")
    projection = reader.admit_arm(tmp_path, plan, rows[1], ["a", "b"])
    row = rows[3]
    if mutation == "no_projection": projection = None
    elif mutation == "copied_points": row["scores"]["FAS"] = .99
    elif mutation == "reuse_pin": row["measurement_scoring_reuse"]["sha256"] = "wrong"
    elif mutation == "resampled": row["metric_files"] = [{"new_FAS_sample": True}]
    elif mutation == "native": projection["products"]["product.tsv"]["sha256"] = "wrong"
    else: projection["candidate"]["membership_constraints"]["sha256"] = "wrong"
    with pytest.raises(ValueError): reader.admit_arm(tmp_path, plan, row, ["a", "b"], projection)


@pytest.mark.parametrize("mutation", ["order", "missing", "source", "digest", "policy", "scores", "graph", "tree", "coverage", "populations", "reuse", "role"])
def test_full_dataset_cannot_ignore_changed_measurement_or_original_roles(tmp_path, monkeypatch, mutation):
    plan, _, cells, _ = admitted_panel(tmp_path, monkeypatch)
    if mutation == "order": cells.reverse()
    elif mutation == "missing": cells.pop()
    elif mutation == "source": cells[4]["source_revision"] = "other"
    elif mutation == "digest": cells[4]["activation"]["candidate_digest"] = "other"
    elif mutation == "policy": cells[4]["activation"]["projection_policy"] = "legacy"
    elif mutation == "scores": cells[3]["scores"][reader.ENDPOINTS[-1]] += .01
    elif mutation == "graph": cells[3]["graph_stages"][0]["sha256"] = "other"
    elif mutation == "tree": cells[3]["rooted_species_tree"]["sha256"] = "other"
    elif mutation == "coverage": cells[3]["coverage"]["pairs"] += 1
    elif mutation == "populations": cells[3]["functional_populations"] = {"changed": True}
    elif mutation == "reuse": cells[3]["assessment_reuse"] = None
    else: cells[2]["role"] = "new_winner"
    with pytest.raises(ValueError): reader.verify_dataset(cells, plan["selection"])


def vector_fixture(passing=None, inactive=()):
    pilot, producer = evidence(passing, inactive)
    selected = reader.planner.selection(pilot, producer)
    cells = [{"label": label, "dataset": dataset, "scores": {key: .7 for key in reader.ENDPOINTS if key.startswith(dataset + ":")}}
        for dataset in reader.DATASETS for label in selected["execution_arms"]]
    return pilot, cells, selected


def test_all_eight_slots_and_distinct_backend_contrasts_are_retained(tmp_path):
    pilot, cells, selected = vector_fixture()
    methods, contrasts = reader.assemble(pilot, cells, selected)
    assert all(len(method["scores"]) == 8 for method in methods)
    assert all(method["full_vector_screen"]["eligible"] for method in methods[1:] if method["role"] == "advancing_candidate")
    assert not methods[3]["full_vector_screen"]["eligible"]
    assert all(value == 0 for value in contrasts["hmmer_versus_builtin"].values())
    assert set(contrasts["hmmer_versus_builtin"]) == set(reader.ENDPOINTS)
    reader.write_table(tmp_path / "scores.tsv", methods, [])
    assert "measurement_control" in (tmp_path / "scores.tsv").read_text()


@pytest.mark.parametrize("passing", [["projection_off"], ["projection_builtin_reselect"], ["projection_hmmer_reselect"]])
def test_stopped_candidates_and_required_controls_never_become_advancing_candidates(passing):
    pilot, cells, selected = vector_fixture(passing)
    methods, _ = reader.assemble(pilot, cells, selected)
    for method in methods:
        label = method["cells"][0]["label"]
        if label not in passing:
            assert not method.get("full_vector_screen", {}).get("eligible", False)
        if label not in selected["execution_arms"]:
            assert len(method["missing_endpoints"]) == 7


def test_applicability_fallback_cannot_relabel_failed_pilot_as_full_vector_pass():
    pilot, cells, selected = vector_fixture([], ["projection_hmmer_reselect"])
    for cell in cells:
        if cell["label"] == "projection_hmmer_reselect": cell["scores"] = {key: .99 for key in cell["scores"]}
    methods, _ = reader.assemble(pilot, cells, selected)
    method = methods[-1]
    assert method["role"] == "applicability_followup"
    assert method["full_vector_screen"]["point_screen_pass"]
    assert not method["full_vector_screen"]["eligible"]


@pytest.mark.parametrize("mutation", ["missing", "order", "unknown", "decline", "support_only"])
def test_original_full_vector_screen_and_complete_cell_inventory_are_preserved(mutation):
    pilot, cells, selected = vector_fixture()
    if mutation == "missing": cells.pop()
    elif mutation == "order": cells.reverse()
    elif mutation == "unknown": cells[1]["scores"]["invented"] = .9
    elif mutation == "decline":
        next(c for c in cells if c["dataset"] == "QfO" and c["label"] == "projection_off")["scores"][reader.ENDPOINTS[5]] = .1
    else:
        control = pilot["cells"][0]["scores"]
        pilot["cells"][1]["scores"] = dict(control)
        next(c for c in cells if c["dataset"] == "QfO" and c["label"] == "projection_off")["scores"][reader.ENDPOINTS[3]] = .9
    if mutation in ("missing", "order", "unknown"):
        with pytest.raises(ValueError): reader.assemble(pilot, cells, selected)
    else:
        methods, _ = reader.assemble(pilot, cells, selected)
        assert not methods[1]["full_vector_screen"]["eligible"]


@pytest.mark.parametrize("revision", [None, "HEAD", "abc123", "f" * 39, "F" * 40])
def test_launch_requires_full_immutable_revision_before_git_lookup(tmp_path, monkeypatch, revision):
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **k: pytest.fail("No mutable revision lookup"))
    with pytest.raises(ValueError, match="immutable"): reader.launch(tmp_path, revision)


def launch_fixture():
    return {"runner_commit": reader.DRIVER, "scientific_revision": followup.SCIENCE, "source": "driver",
        "prepared_inputs": {"path": followup.PREPARED, "bytes": 42, "sha256": "fixed"},
        "independent_reader": {"runner_commit": "a" * 40, "scientific_reader_source": "reader"},
        "datasets": [{"dataset": dataset, "execution_plan": {"path": dataset + "/execution_plan.json", "bytes": 42, "sha256": "fixed"},
            "inference": {"submitted": True, "job_id": str(25000 + i), "planned_cpus": 16, "planned_memory": "128G", "planned_wall_limit": "1-00:00:00"}}
            for i, dataset in enumerate(reader.DATASETS)]}


@pytest.mark.parametrize("mutation", [None, "blob", "science", "driver", "datasets", "unsubmitted", "cpu", "memory", "wall", "job", "duplicate", "prepared", "escape", "reader"])
def test_immutable_launch_binds_two_actual_jobs_resources_paths_and_reader(tmp_path, monkeypatch, mutation):
    entry = launch_fixture()
    if mutation == "science": entry["scientific_revision"] = "other"
    elif mutation == "driver": entry["runner_commit"] = "other"
    elif mutation == "datasets": entry["datasets"].reverse()
    elif mutation == "unsubmitted": entry["datasets"][0]["inference"]["submitted"] = False
    elif mutation == "cpu": entry["datasets"][0]["inference"]["planned_cpus"] = 32
    elif mutation == "memory": entry["datasets"][0]["inference"]["planned_memory"] = "64G"
    elif mutation == "wall": entry["datasets"][0]["inference"]["planned_wall_limit"] = "2-00:00:00"
    elif mutation == "job": entry["datasets"][0]["inference"]["job_id"] = None
    elif mutation == "duplicate": entry["datasets"][1]["inference"]["job_id"] = entry["datasets"][0]["inference"]["job_id"]
    elif mutation == "prepared": entry["prepared_inputs"]["path"] = "other.json"
    elif mutation == "escape": entry["source"] = "../outside"
    elif mutation == "reader": entry["independent_reader"]["runner_commit"] = "HEAD"
    def git(command, **kwargs):
        if command[1] == "cat-file": return "blob" if mutation == "blob" else "commit"
        return json.dumps({"projection_followup_execution": entry}).encode()
    monkeypatch.setattr(reader.subprocess, "check_output", git)
    if mutation:
        with pytest.raises(ValueError): reader.launch(tmp_path, "f" * 40)
    else:
        actual, binding = reader.launch(tmp_path, "f" * 40)
        assert actual == entry and binding["revision"] == "f" * 40


@pytest.mark.parametrize("mutation", [None, "partial", "job", "source", "missing", "role", "failed", "promoted"])
def test_require_complete_ordered_producer_and_original_roles(monkeypatch, mutation):
    selected = chosen()
    plan = {"selection": selected, "driver_revision": reader.DRIVER, "dataset": "QfO"}
    report = {"status": "all_projection_followup_arms_scored", "job_id": "25000", "source_revision": followup.SCIENCE,
        "driver_revision": reader.DRIVER, "dataset": "QfO", "selection": selected,
        "arms": [{"label": label, "dataset": "QfO", "role": selected["roles"][label], "status": followup.SUCCESS} for label in selected["execution_arms"]],
        "defaults_promoted": False, "goal_completed": False}
    entry = {"dataset": "QfO", "inference": {"job_id": "25000"}}
    monkeypatch.setattr(reader.followup, "validate_plan", lambda _: None)
    if mutation == "partial": report["status"] = "running"
    elif mutation == "job": report["job_id"] = "wrong"
    elif mutation == "source": report["source_revision"] = "wrong"
    elif mutation == "missing": report["arms"].pop()
    elif mutation == "role": report["arms"][2]["role"] = "new_winner"
    elif mutation == "failed": report["arms"][2]["status"] = "failed"
    elif mutation == "promoted": report["defaults_promoted"] = True
    if mutation:
        with pytest.raises(ValueError): reader.require_completed(plan, report, entry, selected)
    else: reader.require_completed(plan, report, entry, selected)


def test_saved_functional_numeric_populations_are_distinct_from_eligible_pool(tmp_path):
    pins, details, scores = [], {}, {}
    for metric in ("EC", "GO", "FAS"):
        folder = tmp_path / "qfo_scoring/results" / metric
        folder.mkdir(parents=True)
        raw = folder / (metric + "_fixture_raw.txt.gz")
        with gzip.open(raw, "wt") as stream:
            if metric == "FAS": stream.write("Acc1\tAcc2\tFAS\nA\tB\t0.5\nA\tC\t0.7\n")
            else: stream.write(f"# {metric} Similarities between orthologs from fixture\n# Computing timestamp: today\n# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity\nA\tB\t0.500000\nA\tC\t0.700000\n")
        pins.append(record(raw))
        scores[metric] = .6
        details[metric] = {"assessed_relations": 20 if metric == "FAS" else 2}
        if metric == "FAS":
            write_json(folder / "FAS.json", {"datalink": {"inline_data": {"challenge_participants": [
                {"metric_x": 20, "metric_y": .6, "stderr_y": statistics.stdev([.5, .7]) / 2 ** .5}]}}})
        else: write_json(folder / (metric + ".json"), {"fixture_aggregate": metric, "score": .6})
        pins.append(record(folder / (metric + ".json")))
    populations, checked = reader.functional_populations(tmp_path, scores, details, pins)
    assert populations["FAS"]["sample_pairs"] == 2
    assert populations["FAS"]["reported_eligible_pairs"] == 20
    assert populations["FAS"]["proteins_in_multiple_sample_pairs"] == 1
    assert populations["EC"]["numeric_scored_pairs"] == populations["GO"]["numeric_scored_pairs"] == 2
    assert set(pin["path"] for pin in checked) == set(pin["path"] for pin in pins)
    aggregates_only = [pin for pin in pins if pin["path"].endswith(".json")]
    assert reader.functional_populations(tmp_path, scores, details, aggregates_only)[0] == populations
    with pytest.raises(ValueError, match="absent"): reader.functional_populations(tmp_path, scores, details, [])
    details["GO"]["assessed_relations"] = 3
    with pytest.raises(ValueError, match="count/mean"): reader.functional_populations(tmp_path, scores, details, pins)


def test_existing_reader_output_is_preserved_before_any_launch_lookup(tmp_path):
    with pytest.raises(FileExistsError): reader.run(tmp_path, "f" * 40, tmp_path)


@pytest.mark.parametrize("changed", [False, True])
def test_reader_inventory_rejects_uncommitted_source(tmp_path, monkeypatch, changed):
    (tmp_path / "benchmark_tools").mkdir()
    files = [tmp_path / "benchmark_tools" / name for name in ("fixture.py", reader.PROTOCOL, "hmm_projection_followup_20261019.json")]
    for path in files: path.write_text("frozen\n")
    monkeypatch.setattr(reader, "sources", lambda _: [record(files[0])])
    monkeypatch.setattr(reader.subprocess, "check_output", lambda command, **kwargs: "a" * 40 if command[1] == "rev-parse" else b"changed\n" if changed else b"frozen\n")
    if changed:
        with pytest.raises(ValueError, match="committed"): reader.reporter_inventory(tmp_path)
    else:
        revision, pins = reader.reporter_inventory(tmp_path)
        assert revision == "a" * 40 and len(pins) == 3


def run_fixture(tmp_path, monkeypatch, mutation):
    selected = chosen()
    science, driver, reporter = (tmp_path / name for name in ("science", "driver", "reader"))
    for root in (science, driver, reporter): (root / "benchmark_tools").mkdir(parents=True)
    slurm = science / "benchmark_tools/run_profile_weighting_pilot.slurm"
    slurm.write_text("frozen science\n")
    for path in reader.DRIVER_FILES: (driver / path).write_text("frozen driver\n")
    reader_file = reporter / "benchmark_tools/report_projection_followup.py"
    reader_file.write_text("frozen reader\n")
    protocol = json.loads(Path(reader.__file__).with_name(reader.PROTOCOL).read_text())
    write_json(reporter / "benchmark_tools" / reader.PROTOCOL, protocol)
    prepared_path = tmp_path / followup.PREPARED
    prepared_path.parent.mkdir(parents=True)
    write_json(prepared_path, {"runtime": {"files": []}, "inputs": []})
    prepared_pin = record(prepared_path)
    audit = {"selection": selected, "pilot_readback": record(reader_file), "pilot_plan": record(reader_file),
        "scientific_source": str(science), "runtime": {"files": []}, "tool_path": "/tools",
        "scientific_sources": [record(slurm)], "hmmer_tools": {}}
    entry = {"runner_commit": reader.DRIVER, "scientific_revision": followup.SCIENCE, "source": "driver",
        "selection": selected, "prepared_inputs": {**prepared_pin, "path": followup.PREPARED},
        "independent_reader": {"runner_commit": "a" * 40, "scientific_reader_source": "reader"}, "datasets": []}
    plans, native_callbacks, activation_callbacks, score_callbacks = {}, {}, {}, {}
    for i, dataset in enumerate(reader.DATASETS):
        folder = tmp_path / dataset
        folder.mkdir()
        plan, rows, _ = panel_fixture(folder, monkeypatch, dataset)
        native_callbacks[dataset] = reader.replay_native
        activation_callbacks[dataset] = reader.validate_activation
        score_callbacks[dataset] = reader.readback_scores
        plan.update(**audit, driver_revision=reader.DRIVER, repo=str(tmp_path), driver_root=str(driver),
            driver_sources=[record(driver / path) for path in reader.DRIVER_FILES], prepared_inputs=prepared_pin)
        plan["inputs"].append(prepared_pin)
        if dataset == "QfO":
            relative = "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
            recovery = tmp_path / relative
            recovery.mkdir(parents=True)
            write_json(recovery / "execution_plan.json", json.loads((folder / relative / "execution_plan.json").read_text()))
            original = json.loads((folder / relative / "results.json").read_text())
            original["execution_plan"] = record(recovery / "execution_plan.json")
            write_json(recovery / "results.json", original)
            rows[0]["control_scoring_reuse"] = record(recovery / "results.json")
            plan["inputs"] += [record(recovery / name) for name in ("results.json", "execution_plan.json")]
        for row in rows:
            arm = Path(plan["output"]) / row["label"]
            graph, native = reader.arm_commands(plan, arm, row["label"])
            write_json(arm / "preflight.json", {"graph": graph, "native": native})
            row.update(graph_execution=timing(arm, "graph", graph), native_execution=timing(arm, "native", native))
        path = Path(plan["output"]) / "execution_plan.json"
        write_json(path, plan)
        producer = {"status": "all_projection_followup_arms_scored", "execution_plan": record(path), "job_id": str(25000 + i),
            "source_revision": followup.SCIENCE, "driver_revision": reader.DRIVER, "dataset": dataset,
            "selection": selected, "arms": rows, "defaults_promoted": False, "goal_completed": False}
        if mutation == "partial" and i == 0 or mutation == "second_failed" and i == 1: producer["status"] = "running"
        if mutation == "measurement" and i == 0: producer["arms"][3]["metric_vector"][reader.ENDPOINTS[1]] += .01
        write_json(Path(plan["output"]) / "results.json", producer)
        entry["datasets"].append({"dataset": dataset, "execution_plan": {**record(path), "path": str(path.relative_to(tmp_path))},
            "inference": {"job_id": str(25000 + i)}})
        plans[dataset] = plan
    pilot, _, _ = vector_fixture()
    pilot["uncertainty"] = {"scope": "mocked provenance; no new uncertainty inferred"}
    monkeypatch.setattr(reader, "__file__", str(reader_file))
    monkeypatch.setattr(reader, "launch", lambda *args: (entry, {"revision": "f" * 40}))
    monkeypatch.setattr(reader, "reporter_inventory", lambda *args: ("a" * 40, [record(reader_file)]))
    monkeypatch.setattr(reader.followup, "authenticate_prepared", lambda *args: audit)
    monkeypatch.setattr(reader.planner, "authenticate", lambda *args: (pilot, {}, selected, [], {"fixture": "authenticated"}))
    monkeypatch.setattr(reader.pilot_reader, "retained_context", lambda *args: ([], []))
    monkeypatch.setattr(reader.pilot_reader, "terminal", lambda job: {"job_id": job, "state": "COMPLETED"})
    monkeypatch.setattr(reader, "sequential_science", lambda entries: {"scope": "mocked ordered accounting; separate tests use raw sacct"})
    monkeypatch.setattr(reader.followup, "validate_plan", lambda _: None)
    monkeypatch.setattr(reader, "sources", lambda _: [])
    monkeypatch.setattr(reader, "runtime", lambda: {"files": []} if mutation != "runtime" else {"changed": True})
    monkeypatch.setattr(reader.planner, "selected_case", lambda api, repo, dataset: (plans[dataset]["case"], []))
    monkeypatch.setattr(reader.weighting, "ownership", lambda *args: ["a", "b"])
    monkeypatch.setattr(reader, "replay_native", lambda row, *args: native_callbacks[row["dataset"]](row, *args))
    monkeypatch.setattr(reader, "validate_activation", lambda graph, label, *args, **kwargs: activation_callbacks["QfO"](graph, label, *args, **kwargs))
    monkeypatch.setattr(reader, "readback_scores", lambda repo, dataset, *args: score_callbacks[dataset](repo, dataset, *args))
    monkeypatch.setattr(reader.subprocess, "check_output", lambda command, **kwargs: followup.SCIENCE if kwargs["cwd"] == science else reader.DRIVER)
    if mutation == "reader": entry["independent_reader"]["runner_commit"] = "b" * 40
    elif mutation == "selection":
        entry["selection"] = deepcopy(selected)
        entry["selection"]["roles"]["projection_off"] = "new_winner"
    return tmp_path / "readback"


@pytest.mark.parametrize("mutation", [None, "partial", "second_failed", "runtime", "reader", "selection", "measurement"])
def test_full_reader_lifecycle_admits_only_complete_verified_vectors_and_preserves_failures(tmp_path, monkeypatch, mutation):
    output = run_fixture(tmp_path, monkeypatch, mutation)
    if mutation:
        with pytest.raises(ValueError): reader.run(tmp_path, "f" * 40, output)
        if output.exists():
            report = json.loads((output / "summary.json").read_text())
            assert report["status"] == "failed" and not report["accuracy_result_admitted"]
            if mutation == "second_failed": assert len(report["cells"]) == 5
    else:
        report = reader.run(tmp_path, "f" * 40, output)
        assert report["status"] == "projection_full_vectors_independently_read_back" and report["accuracy_result_admitted"]
        assert len(report["cells"]) == 10 and len(report["datasets"]) == 2
        assert len(report["methods"]) == 5 and all(len(method["scores"]) == 8 for method in report["methods"])
        assert not report["defaults_promoted"] and not report["goal_completed"]
        assert record(output / "scores.tsv") == report["score_table"]
    if output.exists():
        with pytest.raises(FileExistsError): reader.run(tmp_path, "f" * 40, output)


@pytest.mark.parametrize("mutation", [None, "overlap", "reversed", "missing", "failed", "unknown_time"])
def test_scientific_jobs_are_independently_verified_sequential_not_quiet(tmp_path, monkeypatch, mutation):
    entries = [{"inference": {"job_id": "25000"}}, {"inference": {"job_id": "25001"}}]
    rows = ["25000|COMPLETED|0:0|2026-10-19T10:00:00|2026-10-19T11:00:00",
            "25001|COMPLETED|0:0|2026-10-19T11:00:00|2026-10-19T12:00:00"]
    if mutation == "overlap": rows[1] = rows[1].replace("11:00:00", "10:59:59")
    elif mutation == "reversed": entries.reverse()
    elif mutation == "missing": rows.pop()
    elif mutation == "failed": rows[1] = rows[1].replace("COMPLETED|0:0", "FAILED|1:0")
    elif mutation == "unknown_time": rows[1] = rows[1].replace("11:00:00", "Unknown")
    text = "JobIDRaw|State|ExitCode|Start|End\n" + "\n".join(rows) + "\n"
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **k: text)
    if mutation:
        with pytest.raises(ValueError): reader.sequential_science(entries)
    else:
        result = reader.sequential_science(entries)
        assert [row["JobIDRaw"] for row in result["accounting"]] == ["25000", "25001"]
        assert "quiet host" in result["scope"]
