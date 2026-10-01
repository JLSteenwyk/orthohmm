"""Archive and independently reproduce admitted QfO parameter arithmetic offline."""

import argparse
import ast
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tarfile

RUNNER = "benchmark_tools/bundle_qfo_parameter_component.py"
SUMMARY = "benchmark_tools/results/qfo_private_cpm_parameter_result_22395.json"
READBACK = "benchmark_tools/results/qfo_parameter_complete_export_20261001/readback.json"
REPRODUCER = "reproduce_qfo_parameter_uncertainty.py"
REPRODUCER_SHA = "6b337c936bf1cf6896fea06de6e1924092c1dc6e5b3871ad621edbdc2138ddf8"
ANALYSIS_SHA = "608ac7512993280190350d15c77aab56b26770e885dcae0bba62c88354d0a6ca"
FIGURE_NAMES = ("scores.tsv", "scores.md", "intervals.tsv", "qfo_parameter_neighborhood.png",
                "qfo_parameter_neighborhood.pdf", "qfo_parameter_neighborhood.svg")
SOURCES = (REPRODUCER, "run_private_cpm_parameter_uncertainty.py",
           "bootstrap_qfo_parameter_neighborhood.py", "audit_qfo_parameter_swiss.py",
           "bootstrap_qfo_swiss_stages.py", "audit_qfo_swiss_counts.py")
SCOPE = "admitted_qfo_parameter_numerical_component_not_native_reproduction"


def identity(content):
    return dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def safe_name(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or ".." in path.parts or str(path) != name or "\\" in name:
        raise ValueError("Unsafe component member path")
    return name


def direct(path):
    path = Path(path)
    if not path.is_absolute() or path.resolve() != path or path.is_symlink():
        raise ValueError("Require direct absolute component path")
    return path


def equal_identity(content, ref):
    if identity(content) != {key: ref[key] for key in ("bytes", "sha256")}:
        raise ValueError("Component evidence bytes differ")


def git_bytes(repo, revision, name):
    return subprocess.check_output(["git", "-C", str(repo), "show", revision + ":" + safe_name(name)])


def semantic_checks(payloads):
    analysis = json.loads(payloads["data/analysis.json"])
    summary = json.loads(payloads["data/summary.json"])
    reproduction = json.loads(payloads["data/reproduction.json"])
    readback = json.loads(payloads["data/export_readback.json"])
    exported = json.loads(payloads["data/export_manifest.json"])
    if (summary["status"] != "qfo_private_cpm_parameter_result_summary"
            or analysis["status"] != "corrected_qfo_parameter_uncertainty_audited"
            or analysis["source"]["sha256"] != ANALYSIS_SHA
            or any(analysis[key] is not True for key in
                   ("scientific_inputs_admitted", "uncertainty_admitted", "complete_panel", "protocol_controls_match"))
            or analysis["controlled_timing"] is not False or analysis["publication_ready"] is not False
            or summary["controlled_timing"] is not False or summary["publication_ready"] is not False
            or summary["complete_panel"] is not True or summary["endpoint_count"] != 18
            or summary["estimated_contrasts"] != 6 or summary["replicates"] != 100000
            or summary["seed"] != 20260925 or summary["multiplicity_endpoints"] != 18
            or summary["source"] != analysis["source"]
            or summary["point_estimates"] != analysis["point_estimates"]
            or summary["families"] != analysis["families"]
            or summary["comparisons"] != [{k: row[k] for k in ("candidate", "reference", "status", "metrics")}
                                         for row in analysis["comparisons"]]
            or (analysis["replicates"], analysis["seed"], analysis["multiplicity_endpoints"],
                analysis["estimated_contrasts"], analysis["quantile_method"], analysis["alpha"])
                != (100000, 20260925, 18, 6, "linear", .05)):
        raise ValueError("Wrong admitted complete parameter design or summary")
    equal_identity(payloads["data/analysis.json"], summary["input"])
    if (reproduction["status"] != "qfo_parameter_uncertainty_numerically_reproduced"
            or reproduction["input"] != summary["input"] or reproduction["endpoints"] != 18
            or reproduction["planned_endpoints"] != 18 or reproduction["absolute_tolerance"] != 1e-12
            or reproduction["source"]["sha256"] != REPRODUCER_SHA
            or reproduction["publication_ready"] is not False
            or reproduction["numpy_version"] != analysis["numpy_version"]):
        raise ValueError("Wrong independent reproduction binding")
    equal_identity(payloads["sources/" + REPRODUCER], reproduction["source"])
    equal_identity(payloads["data/reproduction.json"], summary["reproduction"])
    if (readback["status"] != "qfo_parameter_export_output_readback"
            or readback["input"] != summary["input"] or readback["reproduction"] != summary["reproduction"]
            or readback["historical_exporter_binding"] is not None or readback["complete_panel"] is not True
            or readback["planned_endpoints"] != 18 or readback["estimated_endpoints"] != 18
            or readback["publication_ready"] is not False
            or exported["status"] != "qfo_parameter_results_exported"
            or exported["input"] != readback["input"] or exported["reproduction"] != readback["reproduction"]
            or exported["outputs"] != readback["outputs"] or exported["complete_panel"] is not True
            or exported["historical_exporter_binding"] is not None or exported["estimated_endpoints"] != 18
            or exported["planned_endpoints"] != 18 or exported["publication_ready"] is not False
            or tuple(Path(ref["path"]).name for ref in exported["outputs"]) != FIGURE_NAMES):
        raise ValueError("Wrong complete export/readback binding")
    equal_identity(payloads["data/export_manifest.json"], readback["manifest"])
    for ref, name in zip(exported["outputs"], FIGURE_NAMES):
        equal_identity(payloads["figures/" + name], ref)
    for name, ref in (("scientific_protocol.md", analysis["protocol"]),
                      ("execution_protocol.md", analysis["execution_protocol"]), ("plan.json", analysis["plan"])):
        equal_identity(payloads["data/" + name], ref)
    for name in SOURCES:
        ref = analysis["source"] if name == "run_private_cpm_parameter_uncertainty.py" else next(
            ref for ref in analysis["helpers"] if Path(ref["path"]).name == name)
        equal_identity(payloads["sources/" + name], ref)
    return analysis


def verify(directory, manifest_sha):
    directory = direct(directory)
    raw = direct(directory / "bundle.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != manifest_sha:
        raise ValueError("Reviewed component manifest changed")
    manifest = json.loads(raw)
    if (manifest["schema_version"] != 1 or manifest["scope"] != SCOPE
            or manifest["publication_ready"] is not False or manifest["native_inference"] is not False
            or not re.fullmatch("[0-9a-f]{40}", manifest["source_commit"])):
        raise ValueError("Wrong component scope")
    payloads = {}
    expected = {"data/analysis.json", "data/summary.json", "data/reproduction.json", "data/export_readback.json",
                "data/export_manifest.json", "data/scientific_protocol.md", "data/execution_protocol.md",
                "data/plan.json", "README.md", "LICENSE.md", "requirements.txt", RUNNER,
                *["sources/" + name for name in SOURCES], *["figures/" + name for name in FIGURE_NAMES]}
    for ref in manifest["files"]:
        name = safe_name(ref["path"])
        if name in payloads or name not in expected:
            raise ValueError("Duplicate or unexpected component member")
        path = direct(directory / name)
        payloads[name] = path.read_bytes()
        equal_identity(payloads[name], ref)
    if set(payloads) != expected:
        raise ValueError("Missing component member")
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
    if actual != expected | {"bundle.json"}:
        raise ValueError("Unexpected component filesystem member")
    analysis = semantic_checks(payloads)
    if (manifest["numpy_version"] != analysis["numpy_version"]
            or payloads["requirements.txt"].decode() != "numpy==" + analysis["numpy_version"] + "\n"):
        raise ValueError("Numerical environment declaration differs")
    return dict(status="qfo_parameter_numerical_component_verified", manifest=identity(raw), files=len(payloads),
        bytes=sum(len(value) for value in payloads.values()), numerical_reproduction_executed=False,
        native_inference=False, publication_ready=False), payloads


def numerical_verifier(source):
    if hashlib.sha256(source).hexdigest() != REPRODUCER_SHA:
        raise ValueError("Unchanged numerical verifier source required")
    # Compile only the existing pure arithmetic entry point and its constants.
    nodes = [node for node in ast.parse(source).body if
             isinstance(node, ast.FunctionDef) and node.name == "verify" or
             isinstance(node, ast.Assign) and len(node.targets) == 1 and
             isinstance(node.targets[0], ast.Name) and node.targets[0].id in {"ARMS", "METRICS"}]
    if len(nodes) != 3:
        raise ValueError("Unexpected numerical entry-point structure")
    import numpy as np
    namespace = {"np": np}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), REPRODUCER, "exec"), namespace)
    return namespace["verify"], np


def reproduce(directory, manifest_sha, output):
    directory, output = direct(directory), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.is_relative_to(directory):
        raise ValueError("Reproduction output must not mutate the component")
    result = dict(status="validation_failed", native_inference=False, scientific_inputs_readmitted=False,
                  publication_ready=False)
    try:
        initial, payloads = verify(directory, manifest_sha)
        function, np = numerical_verifier(payloads["sources/" + REPRODUCER])
        analysis = json.loads(payloads["data/analysis.json"])
        if np.__version__ != analysis["numpy_version"]:
            raise ValueError("Use the declared numerical package version")
        endpoints = function(analysis)
        if endpoints != 18 or verify(directory, manifest_sha)[0] != initial:
            raise ValueError("Incomplete arithmetic or component changed during reproduction")
        result.update(status="qfo_parameter_component_numerically_reproduced", component=initial,
            endpoints=18, absolute_tolerance=1e-12, numpy_version=np.__version__, python=sys.version,
            verifier=identity(payloads["sources/" + REPRODUCER]),
            entry_point="unchanged verify function and ARMS/METRICS constants selected with AST",
            limitations=["Same numerical implementation, generator and quantiles; not an independent statistical engine.",
                "Historical absolute paths are provenance only and are not dereferenced.",
                "No raw predictions/reference recount, native inference, controlled timing or rights clearance."])
    except BaseException as error:
        result.update(error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return result


def build(repo, revision, output):
    repo, output = direct(repo), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", revision + "^{commit}"], text=True).strip()
    submitted_runner = git_bytes(repo, commit, RUNNER)
    if Path(__file__).read_bytes() != submitted_runner:
        raise ValueError("Executing builder differs from requested committed source")
    payloads, entries = {}, {}
    def add(name, content, provenance, expected=None):
        safe_name(name)
        if name in payloads:
            raise ValueError("Duplicate export member")
        if expected is not None:
            equal_identity(content, expected)
        payloads[name] = content
        entries[name] = dict(path=name, **identity(content), provenance=provenance)
    def committed(name, source, expected=None):
        add(name, git_bytes(repo, commit, source), dict(git_commit=commit, git_path=source), expected)
    committed("data/summary.json", SUMMARY)
    committed("data/export_readback.json", READBACK)
    summary, readback = (json.loads(payloads[name]) for name in ("data/summary.json", "data/export_readback.json"))
    for name, ref, required in (("data/analysis.json", summary["input"], "benchmarks/work/qfo_private_cpm_parameter_uncertainty_20261001.json"),
        ("data/export_manifest.json", readback["manifest"], "benchmark_tools/results/qfo_parameter_complete_export_20261001/manifest.json")):
        path = direct(repo / required)
        if ref["path"] != str(path):
            raise ValueError("Unexpected local-only evidence path")
        add(name, path.read_bytes(), dict(historical_record=ref), ref)
    analysis = json.loads(payloads["data/analysis.json"])
    committed("data/reproduction.json", "benchmark_tools/results/qfo_private_cpm_parameter_reproduction_20261001.json", summary["reproduction"])
    for ref, name in zip(readback["outputs"], FIGURE_NAMES):
        committed("figures/" + name, "benchmark_tools/results/qfo_parameter_complete_export_20261001/" + name, ref)
    for name, ref in (("scientific_protocol.md", analysis["protocol"]),
                      ("execution_protocol.md", analysis["execution_protocol"]), ("plan.json", analysis["plan"])):
        source = Path(ref["path"]).relative_to(repo).as_posix()
        committed("data/" + name, source, ref)
    for name in SOURCES:
        committed("sources/" + name, "benchmark_tools/" + name)
    for name in (RUNNER, "LICENSE.md"):
        committed(name, name)
    add("requirements.txt", ("numpy==" + analysis["numpy_version"] + "\n").encode(), dict(derived_from="data/analysis.json"))
    readme = ("# QfO Parameter Numerical Component\n\n"
        "Preserves the complete admitted report, original evidence paths and exported figure/table bytes.\n"
        "Those paths are provenance, not relocated raw-data dependencies.\n\n"
        "With the externally reviewed bundle.json SHA256, from outside the checkout:\n\n"
        "```sh\npython -I -B benchmark_tools/bundle_qfo_parameter_component.py verify "
        "--directory /absolute/component --manifest-sha256 REVIEWED_SHA256\n"
        "python -I -B benchmark_tools/bundle_qfo_parameter_component.py reproduce "
        "--directory /absolute/component --manifest-sha256 REVIEWED_SHA256 --output /absolute/new-result.json\n```\n\n"
        "The reproduction requires the exact NumPy version in requirements.txt. No package is installed automatically.\n"
        "Only the original verifier's pure arithmetic function and constants are loaded using AST; no project harness imports.\n"
        "All 18 planned endpoints and raw family confusion counts are retained.\n"
        "Not native inference/reference recount, hermetic dependency restoration, controlled timing, independent biological "
        "validation, third-party rights clearance or a complete publication release.\n")
    add("README.md", readme.encode(), dict(generated_by=RUNNER))
    semantic_checks(payloads)
    if Path(__file__).read_bytes() != submitted_runner:
        raise ValueError("Builder source changed during export preparation")
    manifest = dict(schema_version=1, scope=SCOPE, source_commit=commit, native_inference=False,
        publication_ready=False, numpy_version=analysis["numpy_version"], files=[entries[name] for name in sorted(entries)])
    output.mkdir(parents=True, exist_ok=False)
    for name, content in payloads.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(content)
    with (output / "bundle.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    sha = identity((output / "bundle.json").read_bytes())["sha256"]
    return verify(output, sha)[0]


def archive(directory, manifest_sha, output):
    directory, output = direct(directory), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.is_relative_to(directory):
        raise ValueError("Archive must not mutate component")
    checked, payloads = verify(directory, manifest_sha)
    payloads["bundle.json"] = (directory / "bundle.json").read_bytes()
    with output.open("xb") as stream, gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode="w") as handle:
            for name, content in sorted(payloads.items()):
                info = tarfile.TarInfo(name)
                info.size, info.mode, info.mtime = len(content), 0o644, 0
                handle.addfile(info, io.BytesIO(content))
    with tarfile.open(output, "r:gz") as handle:
        members = handle.getmembers()
        if len(members) != len(payloads) or {m.name for m in members} != set(payloads) or any(not m.isfile() for m in members):
            raise ValueError("Unexpected archive members")
        for member in members:
            if handle.extractfile(member).read() != payloads[member.name]:
                raise ValueError("Archive member differs")
    if verify(directory, manifest_sha)[0] != checked:
        raise ValueError("Component changed during archiving")
    return dict(status="qfo_parameter_component_archive_verified", archive=identity(output.read_bytes()),
        members=len(payloads), component=checked, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--repo", type=Path, required=True)
    builder.add_argument("--revision", required=True)
    builder.add_argument("--output", type=Path, required=True)
    for name in ("verify", "reproduce", "archive"):
        command = commands.add_parser(name)
        command.add_argument("--directory", type=Path, required=True)
        command.add_argument("--manifest-sha256", required=True)
        if name != "verify":
            command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "build":
        result = build(args.repo, args.revision, args.output)
    elif args.command == "verify":
        result = verify(args.directory, args.manifest_sha256)[0]
    else:
        result = globals()[args.command](args.directory, args.manifest_sha256, args.output)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
