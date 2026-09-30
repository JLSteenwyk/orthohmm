"""Admit an unscored recovered seed with explicitly amended refinement runtime."""

import argparse
import array
import ast
from contextlib import ExitStack
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import resource
import struct
import subprocess
import sys
from types import SimpleNamespace

from benchmark_tools.probe_cpm_partition_parser import check, record
from benchmark_tools.probe_cpm_private_runtime import unique_records

NATIVE = "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
PREFLIGHT = "benchmarks/results/qfo_cpm_checkpoint_recovery_preflight_v1"
FAILED = "benchmarks/results/qfo_cpm_checkpoint_recovery_admission_v1"
EXECUTOR = "benchmarks/work/cpm_checkpoint_recovery_v1_20260923"
COMMIT = "0637916c14a81e7a3b31f5aeb5fadfc52b103079"
HISTORICAL_PYTHON = "/home/bizon/anaconda3/bin/python"
CONTRACT = "benchmarks/work/cpm_checkpoint_admission_v1_20260923/benchmark_tools/admit_cpm_checkpoint_recovery.py"
CONTRACT_SHA = "76a49ce9cc1b65a3839fffdb5a9658220e8bf4f07756f7649573173b4d3d6a02"
PROTOCOL = "benchmark_tools/results/QFO_CPM_HELPER_RECOVERY_ADMISSION_PROTOCOL_20260930.md"
PINS = {
    NATIVE + "/status.json": "a7c0e4765629c76b9462343c4bbce0880f0831ff239468b11c521d774ad2d1fd",
    PREFLIGHT + "/status.json": "bde225bce51355b27396afe74bcc9c25db4b6169d0e9755da4284a575c83e00e",
    FAILED + "/status.json": "ddaf1344fbf1585c6fc2e0ee6351174afafb44e6950ad7c851af75b7812feb10",
    "benchmark_tools/results/qfo_cpm_helper_refinement_readback_20260930.json":
        "dc739ca221f2095a746dffe1887ca9a07800eaebd2e0f3bdf577a9173d3e8c02",
    "benchmark_tools/results/qfo_cpm_helper_environment_readback_20260930.json":
        "8b8fbe5ebc1426956ad1e946aac0b243b8cc4617f714d9e489c94238ac810b91",
}
GRAPH = dict(vertices=984137, edges=25501180, directed=False,
    ordered_endpoints_sha256="f182b9e9b23f9158e1c36525bf546b58a4a0cf1f67928d3b829ab822833a0382",
    ordered_weights_sha256="6c978e4b078c167c33eb9c69fbedef4e239475e2d93c517a296a964806bcb042")
CONSTRUCTOR_SHA = "93b36aad4916b85195c4dcc7dccb37d944cb911f7b051cf989ad1479729fded4"


def load(item):
    check([item])
    return json.loads(Path(item["path"]).read_bytes())


def original_contracts(source):
    if record(source)["sha256"] != CONTRACT_SHA:
        raise ValueError("Changed original admission contract")
    nodes = [node for node in ast.parse(source.read_text()).body
             if isinstance(node, ast.FunctionDef) and node.name in ("completed", "parent_contract")]
    if len(nodes) != 2 or any(node.decorator_list for node in nodes):
        raise ValueError("Missing plain-data historical contracts")
    # Reuse exact historical checks without importing their scientific dependencies.
    namespace = dict(csv=csv, io=io, math=math, record=record,
                     sys=SimpleNamespace(executable=HISTORICAL_PYTHON))
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"), namespace)
    return namespace["completed"], namespace["parent_contract"]


def scheduler_failure(accounting, job, cpus):
    rows = [row for row in csv.DictReader(io.StringIO(accounting), delimiter="|") if row["JobID"] == job]
    if (len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "AllocCPUS", "NodeList"))
            != ("FAILED", "1:0", str(cpus), "bizon")):
        raise ValueError("Historical failure state changed")
    return rows[0]


def npy_header(stream, dtype):
    if stream.read(6) != b"\x93NUMPY":
        raise ValueError("Wrong graph array magic")
    version = stream.read(2)
    if version not in (b"\x01\x00", b"\x02\x00"):
        raise ValueError("Unsupported graph array version")
    length_bytes = stream.read(2 if version == b"\x01\x00" else 4)
    if len(length_bytes) != (2 if version == b"\x01\x00" else 4):
        raise ValueError("Truncated graph header length")
    length = struct.unpack("<H" if version == b"\x01\x00" else "<I", length_bytes)[0]
    if length > 65536:
        raise ValueError("Unbounded graph header")
    raw = stream.read(length)
    if len(raw) != length:
        raise ValueError("Truncated graph header")
    header = ast.literal_eval(raw.decode("latin1"))
    if (not isinstance(header, dict) or set(header) != {"descr", "fortran_order", "shape"} or header["descr"] != dtype
            or header["fortran_order"] is not False or not isinstance(header["shape"], tuple)
            or len(header["shape"]) != 1 or type(header["shape"][0]) is not int or header["shape"][0] < 0):
        raise ValueError("Unexpected typed graph array")
    return header["shape"][0]


def stream_graph(payload, vertices, chunk_size=100000):
    if chunk_size < 1 or vertices < 1:
        raise ValueError("Invalid graph bounds")
    endpoints, weights_hash, constructor = [hashlib.sha256() for _ in range(3)]
    with ExitStack() as stack:
        streams = [stack.enter_context((payload / (name + ".npy")).open("rb"))
                   for name in ("sources", "targets", "weights")]
        lengths = [npy_header(stream, dtype) for stream, dtype in zip(streams, ("<i4", "<i4", "<f8"))]
        if len(set(lengths)) != 1:
            raise ValueError("Graph array lengths differ")
        for start in range(0, lengths[0], chunk_size):
            count = min(chunk_size, lengths[0] - start)
            blocks = [stream.read(count * width) for stream, width in zip(streams, (4, 4, 8))]
            if any(len(block) != count * width for block, width in zip(blocks, (4, 4, 8))):
                raise ValueError("Truncated graph array")
            values = []
            for block, code in zip(blocks, ("i", "i", "d")):
                typed = array.array(code)
                if typed.itemsize != (8 if code == "d" else 4):
                    raise ValueError("Unsupported native array width")
                typed.frombytes(block)
                if sys.byteorder != "little":
                    typed.byteswap()
                values.append(typed)
            canonical_pairs, raw_pairs = bytearray(count * 16), bytearray(count * 8)
            for i, (left, right, weight) in enumerate(zip(*values)):
                if not (0 <= left < vertices and 0 <= right < vertices) or not math.isfinite(weight):
                    raise ValueError("Invalid graph endpoint or weight")
                struct.pack_into("<qq", canonical_pairs, i * 16, min(left, right), max(left, right))
                struct.pack_into("<ii", raw_pairs, i * 8, left, right)
            endpoints.update(canonical_pairs)
            constructor.update(raw_pairs)
            weights_hash.update(blocks[2])
        if any(stream.read(1) for stream in streams):
            raise ValueError("Trailing graph array bytes")
    return dict(vertices=vertices, edges=lengths[0], directed=False,
                ordered_endpoints_sha256=endpoints.hexdigest(), ordered_weights_sha256=weights_hash.hexdigest()), constructor.hexdigest()


def optimizer_contract(snapshot, boundary, adapter, parity, graph, payload, launcher):
    metadata = dict(cpm_resolution=.12, seed=4, include_isolates=True, output_directory=str(payload.parent))
    overrides = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                     OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if (snapshot["status"] != "before_native_clustering" or snapshot["accuracy_evaluated"] is not False
            or snapshot["metadata"] != metadata or snapshot["cwd"] != str(launcher)
            or len(snapshot["cpu_affinity"]) != 1 or not snapshot["native_libraries"]
            or any(snapshot["environment"][k] != v for k, v in overrides.items())
            or snapshot["inputs"] != [record(payload / name) for name in
                                      ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")]
            or snapshot["python"] != record(HISTORICAL_PYTHON)):
        raise ValueError("Original optimizer context differs")
    for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers"):
        if snapshot["modules"][name] != record(launcher / (name.replace(".", "/") + ".py")):
            raise ValueError("Original optimizer scientific imports differ")
    arguments = dict(initial_membership=None, weights="weight", n_iterations=2, max_comm_size=0,
                     seed=4, kwargs={"resolution_parameter": .12}, partition_type="leidenalg.VertexPartition.CPMVertexPartition")
    if boundary != dict(accuracy_evaluated=False, calls=[dict(arguments=arguments, before=graph,
                       saved=graph, after=graph, status="optimizer_returned")]):
        raise ValueError("Missing exact completed optimizer call")
    expected_adapter = dict(format="python_pairs", accuracy_evaluated=False, calls=[dict(
        status="constructor_returned", shape=[graph["edges"], 2], dtype="int32", n=graph["vertices"],
        directed=False, ordered_input_bytes_sha256=CONSTRUCTOR_SHA)])
    if adapter != expected_adapter:
        raise ValueError("Original constructor input differs")
    if (parity["status"] != "constructor_parity_verified_before_optimizer" or parity["accuracy_evaluated"] is not False
            or parity["fingerprint"] != graph or parity["saved"] != graph
            or any(parity["differences"][k] != dict(different_edges=0, examples=[]) for k in
                   ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor"))
            or parity["differences"]["constructor_dtype"] != "int32"
            or parity["differences"]["constructor_c_contiguous"] is not True):
        raise ValueError("Original exhaustive constructor parity differs")


def partition_coverage(path, index, expected_groups):
    seen, groups, members = bytearray(len(index)), 0, 0
    with path.open() as stream:
        for line in stream:
            tokens = line.split()
            if not tokens:
                continue
            groups += 1
            for token in tokens:
                position = index.get(token)
                if position is None or seen[position]:
                    raise ValueError("Unknown or duplicate partition gene")
                seen[position], members = 1, members + 1
    if groups != expected_groups or members != len(index) or not all(seen):
        raise ValueError("Incomplete partition or wrong group count")
    return dict(genes=members, groups=groups)


def _admit(root, output, protocol_sha):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (900, 900))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    pins = [record(root / name) for name in PINS]
    if [row["sha256"] for row in pins] != list(PINS.values()):
        raise ValueError("Changed admission prerequisites")
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed recovery admission amendment")
    parent, preflight, failed, readback, preparation_readback = [load(item) for item in pins]
    contract = record(root / CONTRACT)
    completed, parent_contract = original_contracts(Path(contract["path"]))
    accounting = subprocess.check_output(["sacct", "-j", "22081,22154,22155", "-X", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    scheduler = completed(accounting, "22154")
    original_failure = scheduler_failure(accounting, "22081_1", 32)
    failed_admission = scheduler_failure(accounting, "22155", 2)
    if (failed["status"] != "recovery_admission_failed" or failed["seed_admitted"] is not False
            or failed["accuracy_evaluated"] is not False or failed["publication_ready"] is not False):
        raise ValueError("Historical failed admission differs")
    executor = root / EXECUTOR
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Original recovery executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/run_cpm_checkpoint_recovery.py")
    stages = parent_contract(parent, "22154", COMMIT, source, root, root / NATIVE)
    if (parent["preflight"] != pins[1] or preflight["status"] != "cpm_checkpoint_preflight_verified_unscored"
            or preflight["preflight_passed"] is not True or parent["runtime_after"] != preflight["runtime"]
            or parent["optimizer"]["saved_graph"] != GRAPH or preflight["saved_graph"] != GRAPH
            or any(original_failure[k] != v for k, v in preflight["original_scheduler"].items())):
        raise ValueError("Original preflight/runtime/graph lineage differs")
    new = load(readback["report"])
    prepared = load(readback["helper_environment"])
    if (readback["status"] != "helper_runtime_refinement_completed_independently_read_back_not_admitted"
            or preparation_readback["status"] != "helper_environment_preparation_independently_read_back"
            or readback["helper_environment"] != preparation_readback["preparation"]
            or new["helper_environment"] != readback["helper_environment"]
            or new["status"] != "private_runtime_refinement_control_completed_not_admitted"
            or type(new["returncode"]) is not int or new["returncode"] != 0 or new["timed_out"] is not False
            or new["refinement_attempts"] != 1 or prepared["refinement_attempts"] != 0
            or prepared["status"] != "helper_complete_private_environment_prepared"):
        raise ValueError("Require distinct completed helper-runtime refinement")
    if any(report[key] is not False for report in (readback, preparation_readback, new, prepared)
           for key in ("seed_admitted", "accuracy_evaluated", "publication_ready")):
        raise ValueError("Unexpected prior scientific admission")
    records = [*pins, protocol, contract, source, record(__file__),
        record(root / "benchmark_tools/probe_cpm_partition_parser.py"),
        record(root / "benchmark_tools/probe_cpm_private_runtime.py"),
        *parent["checked_records"], *preflight["checked_records"], *parent["optimizer"]["checked_records"],
        readback["report"], readback["helper_environment"], *new["checked_records"],
        *new["logs"], *prepared["checked_records"], *readback["source_git_bindings"],
        readback["partition"], readback["child_report"], readback["gene_names"],
        *[stage["output"] for stage in stages], *parent["refinement_reports"],
        *[phase["log"] for phase in parent["phases"]]]
    # Git bindings carry an extra field; evidence identities contain only file fields.
    records = unique_records([{k: row[k] for k in ("path", "bytes", "sha256")} for row in records])
    check(records)
    output.mkdir(parents=True)
    report = dict(status="helper_recovery_admission_running", seed_admitted=False,
        accuracy_evaluated=False, downstream_admitted=False, publication_ready=False,
        checked_records=records, native_attempts=0, accounting=accounting, scheduler=scheduler,
        original_failure=original_failure, historical_failed_admission=failed_admission,
        source=record(__file__), source_report=pins[0], refinement_readback=pins[3], preparation_readback=pins[4])
    with (output / "started.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    try:
        names = (root / NATIVE / "payload/gene_names.txt").read_text().splitlines()
        index = {name: i for i, name in enumerate(names)}
        if len(names) != len(index) or len(index) != GRAPH["vertices"]:
            raise ValueError("Wrong complete gene universe")
        payload, launcher = root / NATIVE / "payload", Path(preflight["context"]["cwd"])
        expected_command = [str(Path(prepared["prefix"]) / "bin/python"), "-B", source["path"],
                            "--root", str(root), "--output", str(Path(readback["child_report"]["path"]).parent),
                            "--mode", "repeat-refinement"]
        expected_environment = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", PYTHONNOUSERSITE="1",
            PYTHONMALLOC="debug", PYTHONFAULTHANDLER="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1")
        if (new["command"] != expected_command or new["cwd"] != str(launcher)
                or new["environment_overrides"] != expected_environment
                or new["runtime_probe"]["gc_enabled"] is not True
                or new["runtime_probe"]["gc_thresholds"] != [700, 10, 10]):
            raise ValueError("Completed refinement command/environment differs")
        graph, constructor = stream_graph(payload, len(index))
        if graph != GRAPH or constructor != CONSTRUCTOR_SHA:
            raise ValueError("Actual saved graph/constructor fingerprint differs")
        evidence = [json.loads((payload / name).read_bytes()) for name in
                    ("worker_before.json", "native_boundary.json", "constructor_adapter.json", "constructor_parity.json")]
        optimizer_contract(*evidence, graph, payload, launcher)
        if parent["optimizer"]["partition"] != record(root / NATIVE / "orthohmm_working_res/orthohmm_edges_clustered.txt"):
            raise ValueError("Recovered optimizer partition differs")
        if any(parent["optimizer"]["partition"][key] != stages[2]["output"][key] for key in ("bytes", "sha256")):
            raise ValueError("Recovered seed is not original optimizer output")
        children = [load(item) for item in parent["refinement_reports"]]
        completed_child = load(readback["child_report"])
        if completed_child != new["result"]:
            raise ValueError("New-runtime completion JSON differs from controller")
        if len(children) != 2:
            raise ValueError("Missing original independent refinement reports")
        for child, filename in zip(children, ("orthogroups_profiles_refined.txt", "refinement_repeat.txt")):
            if (child["output"] != record(root / NATIVE / filename)
                    or any(child["output"][key] != stages[3]["output"][key] for key in ("bytes", "sha256"))
                    or {k: v for k, v in child.items() if k != "output"} !=
                       {k: v for k, v in new["result"].items() if k != "output"}):
                raise ValueError("Original/new-runtime refinement metadata differs")
        if (new["result"]["output"] != readback["partition"] or new["result"]["groups"] != 390845
                or new["result"]["genes"] != len(index) or new["result"]["refinement_directed_hits"] != 0
                or new["result"]["numeric_checkpoint"]["summary"]["species"] != 78
                or new["result"]["numeric_checkpoint"]["summary"]["genes"] != len(index)
                or any(readback["partition"][key] != stages[3]["output"][key] for key in ("bytes", "sha256"))):
            raise ValueError("Recovered refined seed differs from completed new-runtime execution")
        coverage = [dict(label=stage["label"], output=stage["output"],
                        **partition_coverage(Path(stage["output"]["path"]), index, groups))
                    for stage, groups in zip(stages, (316603, 393142, 314274, 390845))]
        repeated = partition_coverage(Path(readback["partition"]["path"]), index, 390845)
        if parent["refinement_comparison"] != dict(genes=len(index), groups=390845, partition_equal=True):
            raise ValueError("Original recorded refinement comparison differs")
        check(records)
        report.update(status="cpm_helper_runtime_recovered_seed_admitted_unscored", seed_admitted=True,
            seed_partition=stages[-1]["output"], stages=stages, coverage=coverage,
            independent_refinement_coverage=repeated, saved_graph=graph, constructor_bytes_sha256=constructor,
            refinement_runtime_amendment=dict(scope="independent refinement only", preparation=readback["helper_environment"],
                completed_control=readback["report"], original_optimizer_runtime_unchanged=True,
                downstream_runtime_authorized=False),
            missing_original_statistics=parent["missing_original_statistics"],
            limitations=["Original failures remain failed; new-runtime independent refinement is an explicit amendment.",
                "Recovered seed only; candidate handoff, phylogeny, conversion and scoring require separate admissions.",
                "No new optimizer/refinement execution, endpoint/default change or claim of repaired corruption.",
                "Unchanged source/data and matching outputs do not establish universal runtime equivalence or safety.",
                "Shared-host recovery history is not controlled end-to-end resource evidence."])
    except BaseException as error:
        report.update(status="helper_recovery_admission_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "status.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


def admit(root, output, protocol_sha):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    try:
        return _admit(root, output, protocol_sha)
    except BaseException as error:
        output.mkdir(parents=True, exist_ok=True)
        if not (output / "status.json").exists():
            with (output / "status.json").open("x") as stream:
                json.dump(dict(status="helper_recovery_preflight_failed", native_attempts=0,
                    error_type=type(error).__name__, error=str(error), seed_admitted=False,
                    accuracy_evaluated=False, downstream_admitted=False, publication_ready=False),
                    stream, indent=2, sort_keys=True)
                stream.write("\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.absolute(), args.protocol_sha256)
