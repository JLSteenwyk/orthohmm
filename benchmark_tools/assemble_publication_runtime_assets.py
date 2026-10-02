"""Assemble pinned inference/reader/tool assets without old workstation trees.

The new externally anchored format is explicitly selected by the integrated
executor. Legacy manifests, scientific settings and admitted runs stay intact.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools import bundle_publication_source as source
from benchmark_tools import acquire_publication_wheels as wheels
from benchmark_tools import prepare_publication_phylogeny_tools as tools
from benchmark_tools import export_publication_readers as readers
from benchmark_tools.audit_frozen_overlay_install import scientific_members
from benchmark_tools.run_integrated_publication_workflow import HARNESS, LOCKS, record, save

INDEX = "ASSEMBLY_INDEX.json"
SCHEMA = "publication_runtime_assets_v1"
READER_REVISION = "487d759d29ba2e88eb83346690d7fde047dddfda"
READER_MANIFEST_SHA = "c2a77c7a2635b9193f4321f96db65028a4d5158facd8e52bb4285fe68804f334"
READER_PINS = {
    "benchmark_tools/audit_accuracy_checkpoint.py": dict(bytes=4044, sha256="b71fa9cf1d2cb954206e61179535d7e22fc39946c4de97af6fa45533f7eb9643"),
    "benchmark_tools/audit_historical_profile_ablation.py": dict(bytes=10140, sha256="ffdafda2b55c9580eccc7497881f0e160450af6c9029e3540009bebba2081b31"),
    "benchmark_tools/audit_installed_orthobench.py": dict(bytes=9525, sha256="728d7b9767768cdbd9cd31665fbb44e5e2d2655752e34d5bd0514cf9fc024ee5"),
    "benchmark_tools/audit_phylogeny_events.py": dict(bytes=8605, sha256="7190ec019f98100a241b6fa52640147e97e08eaa6c4d124bf659ca2e143ebabe"),
    "benchmark_tools/audit_phylogeny_hierarchy.py": dict(bytes=5516, sha256="0569763fcbbd3485a21b87727b2eb3fc11b01a4028421e2af6ccaa724ff998aa"),
    "benchmark_tools/audit_phylogeny_sequences.py": dict(bytes=7315, sha256="5345d10e14c16f04e8de2eea1055bbf6351989f8ab3535bbb606599d92aaa6d8"),
    "benchmark_tools/audit_phylogeny_structure.py": dict(bytes=8196, sha256="321770cd145a0a766e7d9816beca42f01be3d24482f02d01dbb178e9e98920aa"),
    "benchmark_tools/audit_publication_pipeline.py": dict(bytes=4082, sha256="f68766d98d831e52696d70691c00f23f8342d3ddca0730d81444f8e2ca849e38"),
    "benchmark_tools/benchmark_production.py": dict(bytes=7954, sha256="47afa3439fea5e0c6c9ac9d62a1c1f819394ba56ec65a55f7192ac9407f918dc"),
    "benchmark_tools/build_publication_runtime.py": dict(bytes=5270, sha256="34af7d2933a7bfd5b46b2b3af3fcf9f21bd597f533a6d9a23a332e4415286ade"),
    "benchmark_tools/candidate_hit_order_policy.py": dict(bytes=4201, sha256="cc70b884e40fe73a3c25ef9ae60a2133508127323a11382a8c610e0d192ec811"),
    "benchmark_tools/compare_installed_ob_search.py": dict(bytes=8570, sha256="baeae4ec82966f2b7c1440ba70bdc0814e44081cebb88d4314bddadb82d02c24"),
    "benchmark_tools/derive_phylogeny_events.py": dict(bytes=5824, sha256="ea27a18ec49fd9a5b01d556587bc68d05aad49db5e454d8b18735dbf6ba2f14c"),
    "benchmark_tools/orthobench_stage_diagnostics.py": dict(bytes=13747, sha256="1ea4a651fef6765104b082d3b0d17d5c8a6263b8da61afbdab006cd85159b72e"),
    "benchmark_tools/prepare_ob_candidate_neighborhood.py": dict(bytes=7448, sha256="03843ed17fa9c44ea1ce5249782cdf009d8c7db9ceb85c62aa85bc2b4156ec29"),
    "benchmark_tools/prepare_orthobench_factorial.py": dict(bytes=10824, sha256="8363f1f5720424ed2aa79ad073eac3b09d726be0f8c63157293dfb20683f8bdd"),
    "benchmark_tools/probe_installed_ob_clustering.py": dict(bytes=7602, sha256="f1f7abc653e57cf8e8b3c1c87aada90ee138da836edcbef8797f6a66358c830e"),
    "benchmark_tools/probe_installed_ob_graph.py": dict(bytes=8153, sha256="4f5fce7b57ab1177d5606b0d3f081d85c6b1b27bf32aa0cc41de7852ac05c5fb"),
    "benchmark_tools/probe_ob_candidate_order_scores.py": dict(bytes=12101, sha256="e23aa1e8de1ae368624a1efd036924bb42d1d634092922e7eab782cddbc02ff9"),
    "benchmark_tools/probe_ob_canonical_candidates.py": dict(bytes=5724, sha256="c6f55261ce34e17681f793afce5bfeebcbc615c6df585e4bce4f62628e48060e"),
    "benchmark_tools/replay_high_sensitivity.py": dict(bytes=16533, sha256="1e220981be64cc75ba772968b9d53fbaea9638da4fff9c48071fec48d33bd6eb"),
    "benchmark_tools/replay_phylogeny.py": dict(bytes=11445, sha256="e4de6f39c66c0f7e7f63229188a76ece44dbdc9b7771ab1c3b1465cec63aa539"),
    "benchmark_tools/run_installed_orthobench.py": dict(bytes=8001, sha256="f855bd27ec0ba7e7106394992d9a5ea008867099af262e3816f63eb6da341792"),
    "benchmark_tools/run_publication_pipeline.py": dict(bytes=6504, sha256="acab3b802b90ea2c4bc4ebc679f5f64431b1672ad85ca38f7ced212965edec87"),
    "benchmark_tools/run_simulation_generation.py": dict(bytes=7039, sha256="464e5c7b6f0d0430eadeba19ca151e5b1a2b8dc046ba372e99cce48a127158a2"),
    "benchmark_tools/run_simulation_methods.py": dict(bytes=12030, sha256="cccf579c221c51e8d03d512264efc3d0a61f4c0a598863d2423bdbd90f7921be"),
    "benchmark_tools/score_orthobench_partition.py": dict(bytes=6229, sha256="0e24c302a9d12ce8c82a4111ff0321c87303469b3c866c317dc2629ad66f9746"),
    "benchmark_tools/trace_ob_initial_edges.py": dict(bytes=8728, sha256="4c5a49d7dde6b274a91edbbc210e42cb61ec2efd8a7dd466693f0ea8612c9b2a"),
    "benchmark_tools/validate_profile_runtime.py": dict(bytes=4014, sha256="94f6dc4b3fe18e6a7ac7e9a95d4cdf67f98c6d9c5697cf08683e002c200f42d8"),
    "benchmark_tools/verify_simulation_histories.py": dict(bytes=4464, sha256="7a50ac5cb159f4f674b8acd5c40a316b11ac798e7e7fd3a539430ccc223166f1"),
    "benchmark_tools/verify_ygob_validation.py": dict(bytes=9370, sha256="a1d4d6e092e37a134a6307aee8527d654f132aab0a35c4bf93204aa0b1e4e7a6"),
}


def regular(path, maximum=100_000_000):
    path = Path(path).absolute()
    if path.is_symlink() or not path.is_file() or path.stat().st_size > maximum:
        raise ValueError("Require a bounded regular asset: " + str(path))
    return record(path)


def pinned(path, expected):
    path = Path(path).absolute()
    if path.resolve() != path or not path.is_file() or path.stat().st_size != expected["bytes"]:
        raise ValueError("Asset path alias or size differs")
    actual = regular(path)
    if any(actual[key] != expected[key] for key in ("bytes", "sha256")):
        raise ValueError("Frozen asset identity differs: " + str(path))
    return actual


def entries(directory):
    result = []
    for path in sorted(directory.rglob("*")):
        name = path.relative_to(directory).as_posix()
        if name == INDEX:
            continue
        if path.is_symlink():
            target = path.readlink()
            if (target.is_absolute() or not name.startswith("assets/mafft/bin/")
                    or not path.resolve().is_relative_to(directory / "assets/mafft")
                    or not path.resolve().is_file()):
                raise ValueError("Assembly link escapes the MAFFT tree or is broken")
            result.append(dict(path=name, kind="symlink", target=target.as_posix()))
        elif path.is_file():
            item = regular(path)
            result.append(dict(path=name, kind="file", mode=path.stat().st_mode & 0o777,
                               bytes=item["bytes"], sha256=item["sha256"]))
        elif not path.is_dir():
            raise ValueError("Nonregular assembled asset")
    return result


def validate(root, digest):
    root = Path(root).absolute()
    if root.resolve() != root or not root.is_dir():
        raise ValueError("Require a canonical assembly root")
    index = regular(root / INDEX, 2 * 1024 ** 2)
    if index["sha256"] != digest:
        raise ValueError("Assembly index differs from external anchor")
    value = json.loads((root / INDEX).read_bytes())
    if value.get("schema") != SCHEMA or value.get("scientific_revision") != source.SCIENTIFIC_REVISION:
        raise ValueError("Wrong assembly schema/scientific revision")
    rows = value["entries"]
    names = [row["path"] for row in rows]
    if not 0 < len(rows) < 1000 or len(set(names)) != len(names):
        raise ValueError("Invalid assembly inventory")
    for row in rows:
        path = Path(row["path"])
        if path.is_absolute() or ".." in path.parts or path.as_posix() != row["path"] or row["path"] == INDEX:
            raise ValueError("Unsafe assembly inventory path")
        if row["kind"] == "file" and (
                type(row["bytes"]) is not int or not 0 <= row["bytes"] <= 100_000_000
                or row["mode"] not in {0o644, 0o755}
                or not re.fullmatch(r"[0-9a-f]{64}", row["sha256"])):
            raise ValueError("Invalid assembly file metadata")
        if row["kind"] not in {"file", "symlink"}:
            raise ValueError("Invalid assembly member kind")
    if entries(root) != rows:
        raise ValueError("Assembled asset bytes/modes/inventory/links differ")
    inventory_path = root / "assets/benchmark_tools/results/integrated_wheel_elf_20260927.json"
    wheels.pinned(inventory_path, wheels.INVENTORY_SHA, 1024 ** 2)
    selected = wheels.wheel_rows(json.loads(inventory_path.read_bytes()))
    for role, names in wheels.ROLES.items():
        directory = root / ("assets/wheels" if role == "inference" else "reader_wheels")
        if {path.name for path in directory.iterdir()} != {selected[name]["filename"] for name in names}:
            raise ValueError("Assembly wheel role inventory differs")
        for name in names:
            pinned(directory / selected[name]["filename"], selected[name]["expected"])
    tools.inspect_helpers(root / "assets/mafft")
    for row in tools.acquisition.artifacts()[1:]:
        name = Path(row["relative"]).name
        path = root / ("assets/FastTree" if name == "FastTree" else "assets/notices/fasttree/" + name)
        pinned(path, row)
    for role, path in {"inference": root / "assets/benchmark_tools/results/publication_recovery_requirements_20260926.txt",
                       "reader": root / "reader_requirements.txt"}.items():
        if regular(path, 8192)["sha256"] != LOCKS[role]:
            raise ValueError("Assembly historical lock differs")
    for name, digest in HARNESS.items():
        if regular(root / "assets/benchmark_tools" / name)["sha256"] != digest:
            raise ValueError("Assembly inference harness differs")
    if regular(root / "readers/manifest.json", 65536)["sha256"] != READER_MANIFEST_SHA:
        raise ValueError("Assembly frozen reader manifest differs")
    readers.verify(root / "readers")
    if json.loads((root / "readers/manifest.json").read_bytes())["files"] != READER_PINS:
        raise ValueError("Assembly frozen reader catalog differs")
    return dict(status="publication_runtime_assets_verified", manifest=index,
                files=sum(row["kind"] == "file" for row in rows),
                symlinks=sum(row["kind"] == "symlink" for row in rows),
                payload_bytes=sum(row.get("bytes", 0) for row in rows), scientific_revision=value["scientific_revision"],
                publication_ready=False, controlled_timing=False, redistribution_clearance=False)


def validate_for_executor(assets, reader_tree, digest):
    assets, reader_tree = Path(assets).absolute(), Path(reader_tree).absolute()
    if assets.name != "assets" or assets.resolve() != assets or reader_tree != assets.parent / "readers":
        raise ValueError("Require assembly asset/reader paths without aliases")
    result = validate(assets.parent, digest)
    tools.cpu_compatible()
    return result


def preflight(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.resolve() != output:
        raise ValueError("Require canonical fresh assembly output")
    component, prepared, wheelhouse = [path.absolute() for path in
                                       (args.component, args.prepared_tools, args.wheels)]
    for root in (component, prepared, wheelhouse):
        if root.resolve() != root or not root.is_dir() or output.is_relative_to(root):
            raise ValueError("Require canonical input roots and external fresh output")
    verified = source.verify(component, args.manifest_sha256)
    index = json.loads((component / "SOURCE_INDEX.json").read_bytes())
    if index.get("profile") not in {"native-wheels", "native-build"}:
        raise ValueError("Require source profile with frozen wheel/runtime support")
    workflow = component / "workflow"
    inventory_path = workflow / "benchmark_tools/results/integrated_wheel_elf_20260927.json"
    wheels.pinned(inventory_path, wheels.INVENTORY_SHA, 1024 ** 2)
    selected = wheels.wheel_rows(json.loads(inventory_path.read_bytes()))
    inputs = [regular(component / "SOURCE_INDEX.json"), regular(inventory_path), regular(__file__),
              regular(readers.__file__), regular(tools.__file__), regular(wheels.__file__)]
    actual = {}
    for name, row in selected.items():
        path = args.project_wheel.absolute() if name == "orthohmm" else wheelhouse / row["filename"]
        if path.name != row["filename"]:
            raise ValueError("Wrong supplied artifact basename")
        actual[name] = pinned(path, row["expected"])
        inputs.append(actual[name])
    locks = {}
    for role, filename in {"inference": "publication_recovery_requirements_20260926.txt",
                           "reader": "publication_reader_requirements_20260927_v2.txt"}.items():
        path = workflow / "benchmark_tools/results" / filename
        wheels.pinned(path, LOCKS[role], 8192)
        locks[role] = path
        inputs.append(regular(path))
    harness = {}
    for name, digest in HARNESS.items():
        path = workflow / "benchmark_tools" / name
        if regular(path)["sha256"] != digest:
            raise ValueError("Changed frozen inference harness")
        harness[name] = path
        inputs.append(regular(path))
    payloads = readers.closure(lambda name: (workflow / name).read_bytes())
    if {name: readers.identity(data) for name, data in payloads.items()} != READER_PINS:
        raise ValueError("Source component cannot reproduce the frozen reader closure")
    inputs.extend(regular(workflow / name) for name in payloads)
    completion = regular(prepared / "complete.json", 4 * 1024 ** 2)
    if completion["sha256"] != args.prepared_tools_sha256:
        raise ValueError("Prepared tool completion differs from external anchor")
    report = json.loads((prepared / "complete.json").read_bytes())
    if (report["status"] != "private_frozen_phylogeny_tools_prepared" or report["old_installation_required"] is not False
            or report["source_unchanged"] is not True or report["scientific_inference_executed"] is not False
            or len(report["stages"]) != 5 or any(row["returncode"] != 0 for row in report["stages"])
            or report["inventory"] != tools.inventory(prepared / "tools")
            or report["helper_files"] != tools.inspect_helpers(prepared / "tools/mafft")):
        raise ValueError("Prepared tool inventory/admission scope differs")
    inputs.append(completion)
    for row in report["stages"]:
        log = regular(row["log"]["path"])
        if log != row["log"]:
            raise ValueError("Prepared tool stage log differs from anchored receipt")
        inputs.append(log)
    for row in report["inventory"]:
        if row["kind"] == "file":
            inputs.append(regular(row["path"]))
    sources = {row["git_path"]: dict(content=(component / row["path"]).read_bytes(), git_blob=row["git_blob"])
               for row in index["files"] if row["path"].startswith("scientific/")}
    members = scientific_members(Path(actual["orthohmm"]["path"]), sources)
    if len(members) != 33:
        raise ValueError("Require all 33 frozen scientific wheel members")
    return output, component, prepared, verified, actual, locks, harness, payloads, inputs, members


def copy_file(original, target, mode=0o644):
    target.parent.mkdir(parents=True, exist_ok=True)
    expected = regular(original)
    shutil.copyfile(original, target)
    target.chmod(mode)
    actual = regular(target)
    if any(actual[key] != expected[key] for key in ("bytes", "sha256")):
        raise ValueError("Asset changed while copying")


def build(args):
    output, component, prepared, verified, actual, locks, harness, payloads, inputs, members = preflight(args)
    output.mkdir(parents=True)
    save(output / "started.json", dict(inputs=inputs, source=verified, attempts=1,
                                      installation_performed=False, scientific_inference_executed=False))
    bundle = output / "bundle"
    bundle.mkdir()
    try:
        for role, names in wheels.ROLES.items():
            directory = bundle / ("assets/wheels" if role == "inference" else "reader_wheels")
            for name in sorted(names):
                original = Path(actual[name]["path"])
                copy_file(original, directory / original.name)
        copy_file(locks["inference"], bundle / "assets/benchmark_tools/results/publication_recovery_requirements_20260926.txt")
        copy_file(locks["reader"], bundle / "reader_requirements.txt")
        copy_file(component / "workflow/benchmark_tools/results/integrated_wheel_elf_20260927.json",
                  bundle / "assets/benchmark_tools/results/integrated_wheel_elf_20260927.json")
        for name, original in harness.items():
            copy_file(original, bundle / "assets/benchmark_tools" / name)
        readers.write_export(payloads, READER_REVISION, bundle / "readers")
        shutil.copytree(prepared / "tools/mafft", bundle / "assets/mafft", symlinks=True)
        shutil.copytree(prepared / "tools/notices/mafft", bundle / "assets/notices/mafft")
        for path in sorted((prepared / "tools/fasttree").iterdir()):
            target = bundle / ("assets/FastTree" if path.name == "FastTree" else "assets/notices/fasttree/" + path.name)
            copy_file(path, target, mode=0o755 if path.name == "FastTree" else 0o644)
        value = dict(schema=SCHEMA, scientific_revision=source.SCIENTIFIC_REVISION,
                     workflow_revision=verified["workflow_revision"], source_manifest=regular(component / "SOURCE_INDEX.json"),
                     prepared_tool_manifest=regular(prepared / "complete.json"), reader_revision=READER_REVISION,
                     entries=entries(bundle), installation_performed=False, scientific_inference_executed=False,
                     controlled_timing=False, publication_ready=False, redistribution_clearance=False)
        save(bundle / INDEX, value)
        (bundle / INDEX).chmod(0o644)
        digest = regular(bundle / INDEX)["sha256"]
        validation = validate(bundle, digest)
        for item in inputs:
            if regular(item["path"]) != item:
                raise ValueError("Assembly input changed")
        if tools.inventory(prepared / "tools") != json.loads((prepared / "complete.json").read_bytes())["inventory"]:
            raise ValueError("Prepared tool bytes/modes/links changed during assembly")
        source.verify(component, args.manifest_sha256)
        result = dict(status="publication_runtime_assets_assembled", inputs=inputs, source=verified,
                      assembly=validation, scientific_members=members, source_unchanged=True,
                      installation_performed=False, native_code_executed=False, scientific_inference_executed=False,
                      attempts=1, retry=False, controlled_timing=False, publication_ready=False,
                      security_clearance=False, redistribution_clearance=False, limitations=[
                          "Explicitly selected new assembly format; historical asset manifests/admissions remain unchanged.",
                          "The exact project wheel may be source-reconstructed; no original project archive is required.",
                          "Prepared launcher prefix differs; explicit MAFFT_BINARIES selects the assembled helper tree.",
                          "No base installation, tool recompilation, network acquisition or scientific inference here.",
                          "Executor/installed-payload/numerical validation, OS/security/rights and public delivery remain separate.",
                      ])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="runtime_asset_assembly_failed", error_type=type(error).__name__,
                                         error=str(error), inputs=inputs, attempts=1, retry=False, publication_ready=False))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    for name in ("component", "prepared-tools", "wheels", "project-wheel", "output"):
        builder.add_argument("--" + name, type=Path, required=True)
    for name in ("manifest-sha256", "prepared-tools-sha256"):
        builder.add_argument("--" + name, required=True)
    verifier = commands.add_parser("verify")
    verifier.add_argument("root", type=Path)
    verifier.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    result = build(args) if args.command == "build" else validate(args.root, args.manifest_sha256)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
