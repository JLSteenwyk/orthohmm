"""Validate direct native scaling outputs separately from timing and accuracy."""

import json
from copy import deepcopy
from pathlib import Path
import shlex

from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.score_ygob_groups import membership, read_predictions
from benchmark_tools.simulation_method_outputs import unique_path, orthohmm_pairs, orthofinder_pairs
from benchmark_tools.validate_simulation_outputs import NativeOutputFailure, validate_graph_weights
from benchmark_tools.gnu_time_companion import command as time_command, parse as parse_time


def input_universe(dataset):
    owners, species = {}, []
    for item in dataset["inputs"]:
        check(item)
        path = Path(item["path"])
        name = path.stem
        if name in species:
            raise ValueError("Duplicate scaling species")
        species.append(name)
        count = 0
        for seq in SeqIO.parse(path, "fasta"):
            if seq.id in owners or not seq.id or not seq.seq:
                raise ValueError("Invalid scaling sequence universe")
            owners[seq.id] = name
            count += 1
        if not count:
            raise ValueError("Empty scaling proteome")
    if len(owners) != dataset["proteins"] or len(species) != dataset["proteomes"]:
        raise ValueError("Scaling input counts changed")
    return owners, species


def complete_partition(path, format, owners):
    groups = read_predictions(path, format)
    if set(membership(groups)) != set(owners):
        raise NativeOutputFailure("Native partition does not cover the exact input universe")
    return len(groups)


def validate_orthohmm(run, owners, species):
    config, argv = run["configuration"], run["native_argv"]
    metrics_path = Path(config["metrics"])
    metrics = json.loads(metrics_path.read_text())
    expected_command = [argv[0], str(Path(run["cwd"]) / "orthohmm/__main__.py"), *argv[3:]]
    if metrics["command"] != expected_command or metrics["cwd"] != run["cwd"]:
        raise ValueError("Native OrthoHMM command or working directory differs")
    if metrics["status"] != "complete":
        raise NativeOutputFailure("Native OrthoHMM metrics are incomplete")
    counts = metrics["counts"]
    if counts["genes"] != len(owners) or counts["species"] != len(species):
        raise NativeOutputFailure("Native OrthoHMM input counts differ")
    output = Path(config["output"])
    group_path = output / "orthohmm_orthogroups.txt"
    groups = complete_partition(group_path, "named_groups", owners)
    if counts["orthogroups"] != groups:
        raise NativeOutputFailure("Native OrthoHMM group count differs")
    files = [metrics_path, group_path]
    result = {"input_genes": len(owners), "orthogroups": groups}
    if run["native_method"] == "orthohmm_satellite_v2":
        root_path = output / "orthohmm_phylogeny/orthohmm_root_hogs.tsv"
        roots = complete_partition(root_path, "root_hogs", owners)
        pair_path = output / "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
        previous, pairs = None, 0
        for pair in orthohmm_pairs(pair_path, owners):
            if pair[0] >= pair[1] or (previous is not None and pair <= previous):
                raise NativeOutputFailure("Native OrthoHMM pairs are not canonical, unique and sorted")
            previous, pairs = pair, pairs + 1
        if counts["phylogeny_root_hogs"] != roots or counts["phylogeny_ortholog_pairs"] != pairs:
            raise NativeOutputFailure("Native phylogenetic counts differ")
        files += [root_path, pair_path]
        result.update(root_hogs=roots, native_pair_rows=pairs)
    return {**result, "checked_files": [record(p) for p in files]}


def validate_orthofinder(run, owners, species):
    config = run["configuration"]
    output = Path(config["output"])
    log = unique_path(output, "**/Log.txt")
    lines = log.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if "Started OrthoFinder version " in line]
    ends = [i for i, line in enumerate(lines) if line.endswith(" : OrthoFinder run completed")]
    commands = [shlex.split(line.removeprefix("Command Line: ")) for line in lines if line.startswith("Command Line: ")]
    if len(starts) != 1 or not lines[starts[0]].endswith("Started OrthoFinder version 3.1.5"):
        raise ValueError("Unexpected OrthoFinder native version")
    if commands != [run["native_argv"]]:
        raise ValueError("Native OrthoFinder command differs")
    if len(ends) != 1 or ends[0] <= starts[0]:
        raise NativeOutputFailure("Native OrthoFinder completion is missing or ambiguous")
    expected = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in run["dataset"]["inputs"]}
    copied = [record(p) for p in Path(config["copy_inputs_to"]).iterdir() if p.is_file()]
    actual = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in copied}
    if len(copied) != len(expected) or actual != expected:
        raise ValueError("OrthoFinder input copies differ from frozen original basenames/bytes")
    graph = unique_path(output, "**/OrthoFinder_graph.txt")
    validate_graph_weights(graph)
    mapping = unique_path(output, "**/SequenceIDs.txt")
    clusters = unique_path(output, "**/clusters_OrthoFinder_I*.txt_id_pairs.txt")
    groups = read_checkpoint(clusters, mapping, owners)
    orthologues = unique_path(output, "**/Orthologues")
    pairs = orthofinder_pairs(orthologues.parent, owners, species)
    files = [log, graph, mapping, clusters, *sorted(orthologues.glob("Orthologues_*/*__v__*.tsv"))]
    return {"input_genes": len(owners), "checkpoint_groups": len(groups), "native_pair_rows": len(pairs),
            "checked_files": [record(p) for p in files]}


def relocate_evidence(run, roots):
    mapped = deepcopy(run)
    def path(value):
        original = Path(value)
        for source, target in sorted(roots.items(), key=lambda item: len(Path(item[0]).parts), reverse=True):
            if original.is_relative_to(Path(source)):
                return str(Path(target) / original.relative_to(Path(source)))
        raise ValueError("Evidence path has no declared relocation: " + str(value))
    for key in ("output", "metrics", "copy_inputs_to"):
        if key in mapped["configuration"]:
            mapped["configuration"][key] = path(mapped["configuration"][key])
    for row in mapped["dataset"]["inputs"]:
        row["path"] = path(row["path"])
    return mapped, path


def validate(run, measurement, evidence_roots=None):
    companion = run.get("gnu_time")
    expected = run["native_argv"]
    if companion is not None:
        if set(companion) != {"executable", "output"}:
            raise ValueError("Unexpected GNU-time specification")
        expected = time_command(expected, companion["output"], companion["executable"])
    if measurement["command"] != expected or measurement["cwd"] != run["cwd"]:
        raise ValueError("Measured command differs from frozen native command")
    if measurement["exit_code"] != 0 or measurement["timed_out"]:
        raise NativeOutputFailure("Native process did not finish successfully")
    accessible, access_path = relocate_evidence(run, evidence_roots) if evidence_roots else (run, str)
    owners, species = input_universe(accessible["dataset"])
    if run["native_method"] in {"orthohmm_high_sensitivity", "orthohmm_satellite_v2"}:
        result = validate_orthohmm(accessible, owners, species)
    elif run["native_method"] == "orthofinder_full":
        result = validate_orthofinder(accessible, owners, species)
    else:
        raise ValueError("Unknown native scaling method")
    if companion is not None:
        path = Path(access_path(companion["output"]))
        timing = parse_time(path.read_text())
        if timing["exit_status"] != measurement["exit_code"]:
            raise NativeOutputFailure("GNU-time and measured exit statuses disagree")
        result["gnu_time_companion"] = {"source": record(path), "accounting": timing,
                                        "wrapper_in_collector_wall_time": True}
    return {"status": "native_scaling_outputs_checked", "accuracy_evaluated": False,
            "evidence_relocation": {str(k): str(v) for k, v in (evidence_roots or {}).items()},
            "resource_measurements_admitted": False, "source": record(__file__), **result,
            "limitations": ["Caller must independently verify frozen manifest/source/runtime before and after inference.",
                            "Native validity does not establish a quiet host, valid timing, or prediction accuracy.",
                            "OrthoFinder validation expands native relation rows outside the inference timer."]}
