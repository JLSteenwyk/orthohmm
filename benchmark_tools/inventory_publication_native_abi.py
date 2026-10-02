"""Inventory declared ABI requirements of an externally pinned runtime assembly."""

import argparse
import json
from pathlib import Path

from benchmark_tools import assemble_publication_runtime_assets as assembly
from benchmark_tools.inventory_wheel_elf import scan_abi, scan_wheel
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


def inspect(root, manifest_sha256, readelf):
    validated = assembly.validate(root, manifest_sha256)
    root = Path(root)
    readelf = Path(readelf).resolve(strict=True)
    inspector = record(readelf)
    manifest = json.loads((root / assembly.INDEX).read_bytes())
    wheels, tools, unique = [], [], set()
    for entry in manifest["entries"]:
        if entry["kind"] != "file":
            continue
        relative = entry["path"]
        path = root / relative
        if path.suffix == ".whl":
            identity = record(path)
            if identity["sha256"] in unique:
                continue
            unique.add(identity["sha256"])
            wheel = dict(name=path.name, version=None, wheel=identity,
                         notice_candidates=[], native_members=[])
            result = scan_wheel(wheel, readelf, include_abi=True)
            result["assembly_path"] = relative
            wheels.append(result)
        elif relative.startswith("assets/"):
            with path.open("rb") as handle:
                magic = handle.read(4)
            if magic == b"\x7fELF":
                before = record(path)
                tools.append(dict(assembly_path=relative, identity=before, abi=scan_abi(path, readelf)))
                check(before)
    if not wheels or not tools:
        raise ValueError("Require assembled wheel and native-tool evidence")
    check(inspector)
    if assembly.validate(root, manifest_sha256) != validated:
        raise ValueError("Assembly changed during ABI inspection")
    objects = [obj["abi"]["requirements"] for w in wheels for obj in w["objects"]]
    objects += [tool["abi"]["requirements"] for tool in tools]
    versions = sorted({(need["library"], version["name"], version["flags"])
                       for obj in objects for need in obj["version_requirements"] for version in need["versions"]})
    return dict(schema="publication_native_abi_inventory_v1", status="declared_abi_inventory_completed",
                assembly=validated, inspector=inspector, source=record(__file__),
                scanner_source=record(Path(__file__).with_name("inventory_wheel_elf.py")),
                unique_wheels=len(wheels), wheel_elf_objects=sum(len(w["objects"]) for w in wheels),
                native_tool_objects=len(tools), wheels=wheels, tools=tools,
                machines=sorted({obj["header"]["Machine"] for obj in objects}),
                interpreters=sorted({obj["interpreter"] for obj in objects if obj["interpreter"]}),
                required_library_versions=[dict(library=lib, name=name, flags=flags) for lib, name, flags in versions],
                native_code_executed=False, runtime_resolution_verified=False, cross_host_compatibility_verified=False,
                controlled_timing=False, publication_ready=False, security_clearance=False, redistribution_clearance=False,
                limitations=["Declared ELF architecture/interpreter/version needs are necessary constraints, not runtime resolution.",
                    "Wheel tags and a maximum GLIBC version alone do not establish compatibility.",
                    "CPU instruction requirements, unversioned/dlopen/static dependencies and actual host libraries remain separate.",
                    "The base Python/Conda payload and OS are not included in this assembly scan.",
                    "No scientific workload, installation, host contention measurement or dependency upgrade is performed."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--readelf", type=Path, default=Path("/usr/bin/readelf"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    save(args.output, inspect(args.root, args.manifest_sha256, args.readelf))


if __name__ == "__main__":
    main()
