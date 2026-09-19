import json
from pathlib import Path

import igraph
import leidenalg
import numpy as np
import pytest

from benchmark_tools.checked_python_pair_worker import python_pair_constructor
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_leiden_boundary import observe_partition
from benchmark_tools.validate_checked_replay_payload import validate
from benchmark_tools.run_qfo_checked_full_replay import run


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.mark.parametrize("problem", [None, "boundary", "resolution", "input_hash", "metadata", "module", "coverage", "duplicate", "provenance",
                                     "retained", "retained_hash", "retained_path", "retained_coverage"])
@pytest.mark.parametrize("corrected", [False, True, "sequence", "control", "cpm_low", "cpm_high"])
def test_native_result_gate(tmp_path, monkeypatch, problem, corrected):
    root = tmp_path
    executor = root / "executor"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    payload = root / "stage/payload"
    payload.mkdir(parents=True)
    pairs = np.array([[2, 0], [0, 2]], dtype=np.int32)
    for name, data in (("sources", pairs[:, 0]), ("targets", pairs[:, 1]), ("weights", np.ones(2))):
        np.save(payload / (name + ".npy"), data)
    (payload / "gene_names.txt").write_text("a\nb\nc\nd\n")
    output = root / ("replay" if corrected else "output")
    cpm_arm = corrected if corrected in ("control", "cpm_low", "cpm_high") else None
    resolution = {"cpm_low": .08, "cpm_high": .12}.get(cpm_arm, .1)
    metadata = {"cpm_resolution": resolution, "seed": 4, "include_isolates": True, "output_directory": str(output)}
    write_json(payload / "metadata.json", metadata)
    manifest = {"stage": "initial", "index": 0, "output_directory": str(output), "inputs": [record(payload / name) for name in
                ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy", "metadata.json")]}
    write_json(payload.parent / "payload_manifest.json", manifest)
    with python_pair_constructor(igraph, payload / "constructor_adapter.json"):
        graph = igraph.Graph(n=4, edges=pairs, directed=False)
        graph.es["weight"] = [1., 1.]
        with observe_partition(leidenalg, payload):
            leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, weights="weight", seed=4, resolution_parameter=resolution)
    modules = {}
    for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers"):
        path = launcher / (name.replace(".", "/") + ".py")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# fixture\n")
        modules[name] = record(path)
    helpers = []
    for name in ("repeat_qfo_saved_graph.py", "checked_python_pair_worker.py", "probe_leiden_boundary.py", "checked_replay_payload_worker.py"):
        path = executor / "benchmark_tools" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# fixture\n")
        helpers.append(record(path))
    admission = root / "benchmark_tools/results/qfo_checked_repeats_verified_20260917.json"
    write_json(admission, {})
    provenance = {"source": helpers[-1], "helpers": helpers[:3], "manifest": record(payload.parent / "payload_manifest.json"),
                  "inputs": manifest["inputs"], "stage": "initial", "accuracy_evaluated": False, "admission": record(admission)}
    corrected_plan = sequence_plan = None
    if corrected:
        checkpoint = root / "corrected_checkpoint"
        checkpoint.mkdir()
        (checkpoint / "gene_names.txt").write_bytes((payload / "gene_names.txt").read_bytes())
        write_json(checkpoint / "manifest.json", {})
        checkpoint_record = record(checkpoint / "manifest.json")
        names_record = record(checkpoint / "gene_names.txt")
        admission = root / "corrected_admission.json"
        write_json(admission, {"status": "corrected_high_sensitivity_native_evidence_admitted",
            "accuracy_evaluated": False, "checkpoint_manifest": checkpoint_record,
            "content": {"genes": 984137, "species_ownership": {str(i): i for i in range(78)},
                        "checked_records": [names_record]}})
        plan = root / "corrected_plan.json"
        write_json(plan, {"status": "corrected_replay_command_frozen_unrun", "accuracy_evaluated": False,
            "execution_authorized": False, "admission": record(admission), "checkpoint_manifest": checkpoint_record,
            "output_root": str(root), "checked_records": [names_record, checkpoint_record, record(admission)]})
        corrected_plan = record(plan)
        provenance.update(admission=record(admission), corrected_plan=corrected_plan)
        if corrected == "sequence":
            from benchmark_tools.run_sequence_graph_control import graph_command
            variant = dict(checkpoint_manifest=checkpoint_record, gene_names=names_record,
                output_root=str(root), expected_genes=984137, expected_species=78, cpu=32,
                requested_memory_gib=64, cap=None, expected_clustering_calls=["initial", "multipass"],
                expected_stages=["multipass", "multipass_refined"],
                native_command=graph_command(launcher, checkpoint, checkpoint_record["sha256"], root))
            estimated = dict(checkpoint_audit=dict(manifest=checkpoint_record), checkpoint_files=[names_record],
                             estimate=dict(genes=984137, species_slot_extent=78))
            write_json(admission, dict(status="admitted_qfo_graph_payload_estimated", accuracy_evaluated=False,
                graph_launched=False, variants=dict(all_hits=estimated, top100=estimated)))
            source = executor / "benchmark_tools/sequence_graph_evidence.py"
            source.write_text("# fixture\n")
            write_json(plan, dict(status="corrected_sequence_graph_commands_frozen_unrun",
                execution_authorized=False, accuracy_evaluated=False, graph_feasibility_admitted=False,
                variants=dict(all_hits=variant, top100=dict(variant, cap=100)), cwd=str(launcher),
                payload_report=record(admission), source=record(source), helpers=[],
                python=record(__import__("sys").executable),
                checked_records=[names_record, checkpoint_record, record(admission)]))
            sequence_plan, corrected_plan = record(plan), None
            provenance.pop("corrected_plan")
            provenance.update(admission=record(admission), sequence_plan=sequence_plan, sequence_variant="all_hits")
            provenance["helpers"].append(record(source))
    if cpm_arm is not None:
        import importlib
        context_module = importlib.import_module("cpm_replay_context")
        context = {"arm": cpm_arm, "resolution": resolution, "metadata": dict(metadata), "checked_records": []}
        monkeypatch.setattr(context_module, "evidence", lambda *a: context)
        source = executor / "benchmark_tools/cpm_replay_context.py"
        source.write_text("# fixture\n")
        provenance["helpers"].append(record(source))
        provenance["cpm_context"] = context
    worker = {"status": "before_native_clustering", "accuracy_evaluated": False, "metadata": dict(metadata),
              "cwd": str(launcher), "inputs": manifest["inputs"][:4], "cpu_affinity": [0], "modules": modules,
              "native_libraries": [record(igraph._igraph.__file__)], "python": record(__import__("sys").executable), "observer": helpers[0],
              "environment": {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                              "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}
    partition = output / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    partition.parent.mkdir(parents=True)
    partition.write_text("a c\nb\nd\n")
    retained_partition = None
    if problem and problem.startswith("retained"):
        retained = payload.parent / "partition.txt"
        retained.write_bytes(partition.read_bytes())
        if problem == "retained_coverage":
            retained.write_text("a c\nb\n")
        retained_partition = record(retained)
        # A later stage has replaced the live output. Retrospective admission
        # must use this stage's retained artifact, not the last stage's output.
        partition.write_text("later_stage_gene\n")
        if problem == "retained_hash":
            retained.write_text("a b c d\n")
        elif problem == "retained_path":
            retained_partition = record(partition)
    if problem == "boundary":
        boundary = json.loads((payload / "native_boundary.json").read_text())
        boundary["calls"][0]["after"]["ordered_endpoints_sha256"] = "changed"
        write_json(payload / "native_boundary.json", boundary)
    elif problem == "resolution":
        boundary = json.loads((payload / "native_boundary.json").read_text())
        boundary["calls"][0]["arguments"]["kwargs"]["resolution_parameter"] = .09
        write_json(payload / "native_boundary.json", boundary)
    elif problem == "input_hash":
        adapter = json.loads((payload / "constructor_adapter.json").read_text())
        adapter["calls"][0]["ordered_input_bytes_sha256"] = "changed"
        write_json(payload / "constructor_adapter.json", adapter)
    elif problem == "metadata":
        worker["metadata"]["seed"] = 5
    elif problem == "module":
        worker["modules"]["orthohmm.externals"]["sha256"] = "changed"
    elif problem == "coverage":
        partition.write_text("a c\nb\n")
    elif problem == "duplicate":
        partition.write_text("a c\nb\nd a\n")
    elif problem == "provenance":
        provenance["stage"] = "multipass"
    write_json(payload / "checked_payload_provenance.json", provenance)
    write_json(payload / "worker_before.json", worker)
    if problem and problem != "retained":
        with pytest.raises(ValueError):
            validate(payload, manifest, root, executor, corrected_plan=corrected_plan, retained_partition=retained_partition,
                     sequence_plan=sequence_plan, sequence_variant="all_hits" if sequence_plan else None, cpm_arm=cpm_arm)
    else:
        result = validate(payload, manifest, root, executor, corrected_plan=corrected_plan, retained_partition=retained_partition,
                          sequence_plan=sequence_plan, sequence_variant="all_hits" if sequence_plan else None, cpm_arm=cpm_arm)
        assert result["genes"] == 4 and result["groups"] == 3
        if retained_partition is not None:
            assert retained_partition in result["provenance_checked"]


def test_parent_refuses_existing_output(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path)
