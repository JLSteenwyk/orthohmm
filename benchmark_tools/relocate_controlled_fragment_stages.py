"""Package pinned fragment artifacts and replay the unchanged native reader."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.fragment_trace_artifact_access import ArtifactAccess, identity


MAX_BYTES = 64 * 1024 * 1024
SCHEMA = "controlled_fragment_relocated_stages_v1"
CHECKPOINT = "orthofinder_sequence_only"
KERNELS = ("orthofinder_mcl_to_orthogroups", "readback_controlled_fragment_stages",
           "readback_controlled_fragment_stage_table")


def need(condition, message):
    if not condition:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), **identity(path))


def logical(ref):
    return ref.get("absolute_path", ref.get("path"))


def checked(ref):
    path = Path(logical(ref))
    need(identity(path) == {k: ref[k] for k in ("bytes", "sha256")}, "Changed original artifact: " + str(path))
    return path


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("ascii")


def plan(root, report_path, readback_path, readback_sha):
    root = Path(root).resolve(strict=True)
    expected_ref = pin(readback_path)
    need(expected_ref["sha256"] == readback_sha, "Changed expected native readback")
    expected = json.loads(checked(expected_ref).read_text())
    need(expected["status"] == "independent_native_stages_and_na_table_verified"
         and expected["new_inference_or_scoring"] is False and expected["publication_ready"] is False,
         "Invalid expected readback scope")
    need(str(Path(report_path).resolve()) == expected["report"]["path"], "Different stage report")
    report = json.loads(checked(expected["report"]).read_text())
    selection = json.loads(checked(report["selection"]).read_text())
    rows = {}

    def add(ref, role, target=None):
        path = Path(logical(ref))
        need(path.is_absolute() and str(path) == logical(ref) and ".." not in path.parts
             and path.is_relative_to(root) and path.resolve(strict=True).is_relative_to(root),
             "Artifact outside original repository")
        checked(ref)
        row = dict(logical_path=str(path), target=target or str(Path("payloads") / path.relative_to(root)),
                   bytes=ref["bytes"], sha256=ref["sha256"])
        if str(path) in rows:
            previous = rows[str(path)]
            need({k: previous[k] for k in row} == row, "Conflicting artifact bindings")
        else:
            rows[str(path)] = dict(row, roles=[])
        if role not in rows[str(path)]["roles"]:
            rows[str(path)]["roles"].append(role)

    add(expected_ref, "expected_native_readback")
    add(expected["report"], "stage_report")
    add(report["selection"], "selection")
    add(pin(checked(expected["report"]).with_name("stages.tsv")), "stage_table")
    for ref in report["checked_inputs"]:
        add(ref, "stage_input")
    for ref in selection["checked_inputs"]:
        add(ref, "selection_input")

    bindings = {b["seed"]: b["arms"] for b in selection["bindings"]}
    contexts = {(case["seed"], arm) for case in selection["cases"] if case["method"] == CHECKPOINT
                for arm in ("baseline", "fragment")}
    extra_fastas = set()
    for seed, arm in sorted(contexts):
        binding = bindings[seed][arm][CHECKPOINT]
        add(binding["execution"], "checkpoint_execution")
        execution = json.loads(checked(binding["execution"]).read_text())
        run = execution["methods"]["orthofinder_full"]
        need(run["status"] == "process_succeeded" and run["exit_code"] == 0, "Unsuccessful checkpoint parent")
        directory = Path(run["argv"][2])
        refs = [r for r in run["outputs"] if Path(logical(r)).parent == directory
                and Path(logical(r)).suffix == ".fasta"]
        need(refs and {logical(r) for r in refs} == {str(p) for p in directory.glob("*.fasta")},
             "Copied checkpoint FASTAs differ from execution inventory")
        for ref in refs:
            if logical(ref) not in rows:
                extra_fastas.add(logical(ref))
            add(ref, "checkpoint_copied_fasta")

    sources = dict(previous=expected["previous_reader"]["path"], reader=expected["source"]["path"],
                   mcl=str(root / "benchmark_tools" / (KERNELS[0] + ".py")))
    for ref in (expected["previous_reader"], expected["source"], pin(sources["mcl"])):
        add(ref, "unchanged_kernel", "kernels/benchmark_tools/" + Path(ref["path"]).name)
    for name in ("__init__.py", "fragment_trace_artifact_access.py", "relocate_controlled_fragment_stages.py"):
        add(pin(root / "benchmark_tools" / name), "relocation_runtime", "runner/benchmark_tools/" + name)
    values = [dict(row, roles=sorted(row["roles"])) for _, row in sorted(rows.items())]
    total = sum(row["bytes"] for row in values)
    need(total <= MAX_BYTES, "Fragment component exceeds bounded payload size")
    return dict(schema=SCHEMA, artifacts=values, artifact_payload_bytes=total,
                added_checkpoint_fastas=len(extra_fastas), sources=sources, report=expected["report"],
                expected_native_readback=expected_ref, new_inference_or_scoring=False,
                historical_inference_reproduced=False, publication_ready=False,
                path_mapping="Historical logical labels; verified copied files only; no original-path fallback")


def prepare(root, report_path, readback_path, readback_sha, output):
    output = Path(output)
    need(not output.exists() and not output.is_symlink(), "Existing component output")
    manifest = plan(root, report_path, readback_path, readback_sha)
    output.mkdir(parents=True, exist_ok=False)
    for row in manifest["artifacts"]:
        data = checked(dict(path=row["logical_path"], bytes=row["bytes"], sha256=row["sha256"])).read_bytes()
        need(len(data) == row["bytes"] and hashlib.sha256(data).hexdigest() == row["sha256"], "Original changed during copy")
        target = output / row["target"]
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(data)
    ArtifactAccess(output, manifest["artifacts"]).verify()
    with (output / "manifest.json").open("xb") as stream:
        stream.write(json_bytes(manifest))
    return dict(status="pinned_fragment_component_prepared", manifest=pin(output / "manifest.json"),
                artifact_files=len(manifest["artifacts"]), artifact_payload_bytes=manifest["artifact_payload_bytes"],
                added_checkpoint_fastas=manifest["added_checkpoint_fastas"], new_inference_or_scoring=False)


def load_readers(access, sources):
    names = ["benchmark_tools"] + ["benchmark_tools." + name for name in KERNELS]
    missing = object()
    old = {name: sys.modules.get(name, missing) for name in names}
    package = ModuleType("benchmark_tools")
    package.__path__ = []
    sys.modules["benchmark_tools"] = package
    try:
        for name, logical_path in zip(KERNELS, (sources["mcl"], sources["previous"], sources["reader"])):
            physical = access.physical(logical_path)
            spec = importlib.util.spec_from_file_location("benchmark_tools." + name, physical)
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            setattr(package, name, module)
            spec.loader.exec_module(module)
            if name != KERNELS[0]:
                module.Path = access.path
        return getattr(package, KERNELS[2])
    finally:
        for name, original in old.items():
            if original is missing:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = original


def replay(component, manifest_sha):
    component = Path(component).resolve(strict=True)
    manifest_ref = pin(component / "manifest.json")
    need(manifest_ref["sha256"] == manifest_sha, "Changed component manifest")
    manifest = json.loads((component / "manifest.json").read_text())
    need(manifest["schema"] == SCHEMA and all(manifest[k] is False for k in
         ("new_inference_or_scoring", "historical_inference_reproduced", "publication_ready")), "Changed component scope")
    need(sum(row["bytes"] for row in manifest["artifacts"]) == manifest["artifact_payload_bytes"] <= MAX_BYTES,
         "Invalid component payload size")
    access = ArtifactAccess(component, manifest["artifacts"])
    access.verify()
    for name in ("relocate_controlled_fragment_stages.py", "fragment_trace_artifact_access.py"):
        rows = [row for row in manifest["artifacts"] if row["target"] == "runner/benchmark_tools/" + name]
        need(len(rows) == 1 and identity(Path(__file__).with_name(name)) == {k: rows[0][k] for k in ("bytes", "sha256")},
             "Changed relocation runtime")
    reader = load_readers(access, manifest["sources"])
    result = reader.verify(access.path(manifest["report"]["path"]), manifest["report"]["sha256"])
    expected = json.loads(access.path(manifest["expected_native_readback"]["path"]).read_text())
    need(result == expected, "Relocated native readback differs from retained result")
    access.verify()
    import numpy as np
    return dict(status="relocated_native_readback_matches_retained_result", component=str(component), manifest=manifest_ref,
                artifact_files=len(manifest["artifacts"]), artifact_payload_bytes=manifest["artifact_payload_bytes"],
                added_checkpoint_fastas=manifest["added_checkpoint_fastas"], native_readback=result,
                python_version=sys.version.split()[0], numpy_version=np.__version__,
                copied_kernel_sources=[dict(logical_path=p, copied_path=str(access.physical(p)), **identity(access.physical(p)))
                                       for p in manifest["sources"].values()],
                original_path_fallback=False, original_artifacts_modified=False, new_inference_or_scoring=False,
                historical_inference_reproduced=False, publication_ready=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    prepare_parser = modes.add_parser("prepare")
    for name in ("root", "report", "readback", "output"):
        prepare_parser.add_argument("--" + name, type=Path, required=True)
    prepare_parser.add_argument("--readback-sha256", required=True)
    replay_parser = modes.add_parser("replay")
    replay_parser.add_argument("--component", type=Path, required=True)
    replay_parser.add_argument("--manifest-sha256", required=True)
    replay_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), "Existing output")
    if args.mode == "prepare":
        result = prepare(args.root, args.report, args.readback, args.readback_sha256, args.output)
    else:
        result = replay(args.component, args.manifest_sha256)
        with args.output.open("xb") as stream:
            stream.write(json_bytes(result))
    summary = {k: result[k] for k in ("status", "manifest", "artifact_files", "artifact_payload_bytes", "added_checkpoint_fastas")}
    if "native_readback" in result:
        summary.update({k: result["native_readback"][k] for k in ("stage_rows", "selected_cases", "native_contexts")})
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
