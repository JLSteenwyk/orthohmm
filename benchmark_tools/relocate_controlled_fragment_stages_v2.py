"""Complete FASTA ownership bindings; retain the original failed relocation."""

import argparse
import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools import relocate_controlled_fragment_stages as previous
from benchmark_tools.fragment_trace_artifact_access import ArtifactAccess, identity


REVISION = "all_configured_fasta_ownership_v2"


def plan(root, report_path, readback_path, readback_sha):
    root = Path(root).resolve(strict=True)
    manifest = previous.plan(root, report_path, readback_path, readback_sha)
    report = json.loads(previous.checked(manifest["report"]).read_text())
    selection = json.loads(previous.checked(report["selection"]).read_text())
    bindings = {b["seed"]: b["arms"] for b in selection["bindings"]}
    contexts = {(case["method"], case["seed"], arm) for case in selection["cases"]
                for arm in ("baseline", "fragment")}
    rows = {row["logical_path"]: row for row in manifest["artifacts"]}
    added, fasta_reads = set(), 0
    for method, seed, arm in sorted(contexts):
        binding = bindings[seed][arm][method]
        execution = json.loads(previous.checked(binding["execution"]).read_text())
        parent = "orthofinder_full" if method == previous.CHECKPOINT else method
        run = execution["methods"][parent]
        previous.need(run["status"] == "process_succeeded" and run["exit_code"] == 0, "Unsuccessful native context")
        directory = Path(binding["configured"]["argv"][2] if method.startswith("orthohmm_") else
                         binding["configured"]["copy_inputs_from"] if method == "orthofinder_full" else run["argv"][2])
        inventory = run["outputs"] if method == previous.CHECKPOINT else execution["verified_inputs"]["inputs"]
        refs = [ref for ref in inventory if Path(previous.logical(ref)).parent == directory
                and Path(previous.logical(ref)).suffix == ".fasta"]
        names = [previous.logical(ref) for ref in refs]
        previous.need(names and len(set(names)) == len(names) and set(names) == {str(p) for p in directory.glob("*.fasta")},
                      "Configured FASTA ownership differs from retained execution inventory")
        fasta_reads += len(refs)
        for ref in refs:
            path = previous.checked(ref)
            previous.need(path.is_absolute() and ".." not in path.parts and path.is_relative_to(root)
                          and path.resolve(strict=True).is_relative_to(root), "FASTA outside original repository")
            row = dict(logical_path=str(path), target=str(Path("payloads") / path.relative_to(root)),
                       bytes=ref["bytes"], sha256=ref["sha256"])
            if str(path) in rows:
                previous.need({k: rows[str(path)][k] for k in row} == row, "Conflicting FASTA ownership binding")
            else:
                added.add(str(path))
                rows[str(path)] = dict(row, roles=[])
            if "configured_fasta_ownership" not in rows[str(path)]["roles"]:
                rows[str(path)]["roles"].append("configured_fasta_ownership")
    source = previous.pin(root / "benchmark_tools" / Path(__file__).name)
    rows[source["path"]] = dict(logical_path=source["path"], target="runner/benchmark_tools/" + Path(__file__).name,
                               bytes=source["bytes"], sha256=source["sha256"], roles=["relocation_runtime_correction"])
    manifest["artifacts"] = [dict(row, roles=sorted(row["roles"])) for _, row in sorted(rows.items())]
    manifest["artifact_payload_bytes"] = sum(row["bytes"] for row in manifest["artifacts"])
    previous.need(manifest["artifact_payload_bytes"] <= previous.MAX_BYTES, "Fragment component exceeds bounded payload size")
    manifest.update(preparation_revision=REVISION, configured_native_contexts=len(contexts),
                    configured_fasta_reads=fasta_reads,
                    added_configured_fastas=len(added), original_reader_kernels_changed=False)
    return manifest


def prepare(root, report_path, readback_path, readback_sha, output):
    previous.need(not Path(output).exists() and not Path(output).is_symlink(), "Existing component output")
    manifest = plan(root, report_path, readback_path, readback_sha)
    original_plan = previous.plan
    # Reuse the frozen copy writer, changing only its in-memory plan input.
    previous.plan = lambda *args: manifest
    try:
        result = previous.prepare(root, report_path, readback_path, readback_sha, output)
    finally:
        previous.plan = original_plan
    result.update(preparation_revision=REVISION, added_configured_fastas=manifest["added_configured_fastas"])
    return result


def replay(component, manifest_sha):
    component = Path(component).resolve(strict=True)
    manifest = json.loads((component / "manifest.json").read_text())
    previous.need(identity(component / "manifest.json")["sha256"] == manifest_sha
                  and manifest["preparation_revision"] == REVISION and manifest["original_reader_kernels_changed"] is False,
                  "Changed corrected component binding")
    access = ArtifactAccess(component, manifest["artifacts"])
    sources = [r for r in manifest["artifacts"] if r["target"] == "runner/benchmark_tools/" + Path(__file__).name]
    previous.need(len(sources) == 1 and identity(__file__) == {k: sources[0][k] for k in ("bytes", "sha256")},
                  "Changed relocation correction runtime")
    access.physical(sources[0]["logical_path"])
    result = previous.replay(component, manifest_sha)
    result.update(preparation_revision=REVISION, added_configured_fastas=manifest["added_configured_fastas"],
                  original_reader_kernels_changed=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    preparing = modes.add_parser("prepare")
    for name in ("root", "report", "readback", "output"):
        preparing.add_argument("--" + name, type=Path, required=True)
    preparing.add_argument("--readback-sha256", required=True)
    checking = modes.add_parser("replay")
    checking.add_argument("--component", type=Path, required=True)
    checking.add_argument("--manifest-sha256", required=True)
    checking.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    previous.need(not args.output.exists() and not args.output.is_symlink(), "Existing output")
    if args.mode == "prepare":
        result = prepare(args.root, args.report, args.readback, args.readback_sha256, args.output)
    else:
        result = replay(args.component, args.manifest_sha256)
        with args.output.open("xb") as stream:
            stream.write(previous.json_bytes(result))
    summary = {k: result[k] for k in ("status", "manifest", "artifact_files", "artifact_payload_bytes",
                                     "added_checkpoint_fastas", "added_configured_fastas", "preparation_revision")}
    if "native_readback" in result:
        summary.update({k: result["native_readback"][k] for k in ("stage_rows", "selected_cases", "native_contexts")})
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
