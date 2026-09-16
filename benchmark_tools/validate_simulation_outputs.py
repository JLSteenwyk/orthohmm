"""Native completion and provenance gates for frozen simulation predictions."""

import json
import math
from pathlib import Path
import shlex

from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_simulation_generation import verify_file
from benchmark_tools.simulation_method_outputs import unique_path
from benchmark_tools.verify_ygob_validation import verify_records


class NativeOutputFailure(ValueError):
    """Verified process output establishes a scientific/native failure."""


def verify_process(config, record):
    if record["status"] != "process_succeeded" or record.get("exit_code") != 0:
        raise ValueError("Method process did not succeed")
    if record["argv"] != config["argv"]:
        raise ValueError("Inference command differs from frozen manifest")
    root = Path(config["output"]).resolve()
    metrics = Path(config["metrics"]).resolve() if "metrics" in config else None
    observed = set()
    for item in record.get("outputs", []):
        path = Path(item["absolute_path"]).resolve()
        if path in observed or (not path.is_relative_to(root) and path != metrics):
            raise ValueError("Duplicate or out-of-scope inference artifact")
        verify_file(path, item)
        observed.add(path)
    if not observed:
        raise ValueError("Missing inference output inventory")
    actual = {p.resolve() for p in root.rglob("*") if p.is_file()}
    if metrics is not None:
        actual.add(metrics)
    if actual != observed:
        raise ValueError("Inference output file set changed")
    return observed


def input_signature(records):
    result = {}
    for record in records:
        name = Path(record["absolute_path"]).name
        if name in result:
            raise ValueError("Duplicate input basename")
        result[name] = {"path": name, "bytes": record["bytes"], "sha256": record["sha256"]}
    return result


def validate_orthohmm(method, config, record, inputs, manifest):
    inventory = verify_process(config, record)
    metrics_path = Path(config["metrics"])
    metrics = json.loads(metrics_path.read_text())
    harness = metrics["harness"]
    if harness.get("git_commit") != manifest["core_commit"]:
        raise ValueError("OrthoHMM used wrong source revision")
    expected = input_signature(inputs)
    actual = {r["path"]: r for r in harness["input_manifest"]}
    if actual != expected or len(actual) != len(harness["input_manifest"]):
        raise ValueError("OrthoHMM native input inventory mismatch")
    root = Path(manifest["core_root"])
    sources = verify_records(harness["source_manifest"], root)
    expected_sources = {Path(r["absolute_path"]).resolve() for r in manifest["core_sources"]
                        if Path(r["absolute_path"]).suffix == ".py"}
    if sources != expected_sources:
        raise ValueError("OrthoHMM native source inventory mismatch")
    for item in manifest["core_sources"]:
        verify_file(Path(item["absolute_path"]), item)
    if metrics.get("status") != "complete" or harness.get("exit_code") != 0:
        raise NativeOutputFailure("OrthoHMM native completion not confirmed")
    outputs = verify_records(harness["output_manifest"], Path(config["output"]))
    relative = ("orthohmm_orthogroups.txt" if method == "orthohmm_high_sensitivity" else
                "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv")
    native = (Path(config["output"]) / relative).resolve()
    if native not in outputs or not outputs.issubset(inventory):
        raise ValueError("Native predictions absent from verified output manifests")
    return {"completion": "metrics complete and harness exit zero", "native": str(native),
            "metrics": file_record(metrics_path, metrics_path.parent)}


def validate_orthofinder(config, record, inputs):
    inventory = verify_process(config, record)
    root = Path(config["output"])
    log = unique_path(root, "**/Log.txt")
    text = log.read_text()
    lines = text.splitlines()
    starts = [line for line in lines if "Started OrthoFinder version " in line]
    ends = [line for line in lines if line.endswith(" : OrthoFinder run completed")]
    commands = [line.removeprefix("Command Line: ") for line in lines if line.startswith("Command Line: ")]
    if len(starts) != 1 or not starts[0].endswith("Started OrthoFinder version 3.1.5"):
        raise ValueError("OrthoFinder 3.1.5 native completion not confirmed")
    if len(commands) != 1 or shlex.split(commands[0]) != config["argv"]:
        raise ValueError("OrthoFinder native command mismatch")
    copied = list(Path(config["copy_inputs_to"]).glob("*.fasta"))
    if {p.name: file_record(p, p.parent) for p in copied} != input_signature(inputs):
        raise ValueError("OrthoFinder input copies differ")
    if log.resolve() not in inventory:
        raise ValueError("OrthoFinder completion log not inventoried")
    if not ends:
        raise NativeOutputFailure("OrthoFinder 3.1.5 native completion not confirmed")
    if len(ends) != 1 or lines.index(ends[0]) <= lines.index(starts[0]):
        raise ValueError("OrthoFinder completion marker is ambiguous or precedes start")
    graph = unique_path(root, "**/OrthoFinder_graph.txt")
    validate_graph_weights(graph)
    return {"completion": "native version, command and completion marker verified",
            "native_log": file_record(log, root)}


def validate_graph_weights(path):
    """Check numeric entries of the native MCL matrix, not its header/comments."""
    entries = 0
    with path.open() as handle:
        for line in handle:
            if line.lstrip().startswith("#"):
                continue
            for token in line.split():
                if ":" not in token:
                    continue
                index, value = token.split(":", 1)
                if not index.isdigit():
                    raise ValueError("Malformed OrthoFinder graph index")
                if not math.isfinite(float(value)):
                    raise NativeOutputFailure("Nonfinite OrthoFinder graph weight")
                entries += 1
    return entries
