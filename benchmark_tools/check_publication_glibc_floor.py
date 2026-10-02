"""Check a pinned assembly's declared glibc floor, not complete host compatibility."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re


def record(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def version(value):
    if not isinstance(value, str) or not re.fullmatch(r"\d+\.\d+(?:\.\d+)?", value):
        raise ValueError("Require a numeric glibc release")
    numbers = tuple(int(part) for part in value.split("."))
    return numbers + (0,) * (3 - len(numbers))


def evaluate(inventory, host_system, host_machine, host_glibc):
    if (inventory.get("schema") != "publication_native_abi_inventory_v1"
            or inventory.get("status") != "declared_abi_inventory_completed"
            or inventory.get("machines") != ["Advanced Micro Devices X86-64"]):
        raise ValueError("Require the supported completed x86-64 ABI inventory")
    if host_system != "Linux" or host_machine != "x86_64":
        raise ValueError("The assembled assets require Linux x86-64")
    match = re.fullmatch(r"glibc (\d+\.\d+(?:\.\d+)?)", host_glibc or "")
    if not match:
        raise ValueError("Running controller does not report a supported GNU libc release")
    needed = []
    for row in inventory["required_library_versions"]:
        name = row["name"]
        if not isinstance(name, str):
            raise ValueError("Invalid required version name")
        if name.startswith("GLIBC_"):
            version(name[6:])
            needed.append(name[6:])
    if not needed:
        raise ValueError("Missing declared numeric glibc requirements")
    floor = max(needed, key=version)
    if version(match[1]) < version(floor):
        raise ValueError("Host glibc " + match[1] + " is below declared floor " + floor)
    return dict(host_system=host_system, host_machine=host_machine, controller_glibc=host_glibc,
                declared_glibc_floor=floor, glibc_floor_check_passed=True,
                full_host_compatibility_verified=False, runtime_resolution_verified=False,
                reviewed_native_code_executed=False, controlled_timing=False, publication_ready=False,
                limitations=["Numeric GNU libc floor only, including conservatively retained weak requirements.",
                    "Controller confstr reports its libc, not every provider loaded by an asset or private base.",
                    "Loader presence is not loader resolution; non-glibc, CPU, static/dlopen, base/OS and security checks remain.",
                    "This check does not authorize production timing or repeat scientific/native inference."])


def check_host(root, inventory_path, inventory_sha256, assembly_sha256):
    from benchmark_tools.assemble_publication_runtime_assets import validate
    if not isinstance(inventory_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", inventory_sha256):
        raise ValueError("Require an externally retained ABI inventory digest")
    inventory_path = Path(inventory_path).absolute()
    if (inventory_path.is_symlink() or not inventory_path.is_file()
            or inventory_path.stat().st_size > 2 * 1024 ** 2):
        raise ValueError("Require a regular, bounded ABI inventory")
    pinned = record(inventory_path)
    if pinned["sha256"] != inventory_sha256:
        raise ValueError("ABI inventory differs from its external anchor")
    inventory = json.loads(inventory_path.read_bytes())
    validated = validate(root, assembly_sha256)
    recorded = inventory["assembly"]
    if any(recorded.get(key) != validated.get(key) for key in
           ("status", "files", "symlinks", "payload_bytes", "scientific_revision")):
        raise ValueError("ABI inventory describes a different assembly")
    if recorded["manifest"]["sha256"] != validated["manifest"]["sha256"]:
        raise ValueError("ABI inventory describes a different assembly index")
    overrides = ("LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT")
    if any(os.environ.get(key) for key in overrides):
        raise ValueError("Refuse controller loader overrides for glibc observation")
    try:
        glibc = os.confstr("CS_GNU_LIBC_VERSION")
    except (AttributeError, OSError, ValueError) as error:
        raise ValueError("GNU libc observation unavailable") from error
    result = evaluate(inventory, platform.system(), platform.machine(), glibc)
    loaders = []
    for name in inventory["interpreters"]:
        path = Path(name)
        if not path.is_absolute() or not path.is_file() or not os.access(path, os.X_OK):
            raise ValueError("Declared program interpreter is unavailable: " + name)
        loaders.append(dict(declared_path=name, identity=record(path)))
    if not loaders:
        raise ValueError("Missing declared program-interpreter evidence")
    if record(inventory_path) != pinned or validate(root, assembly_sha256) != validated:
        raise ValueError("ABI inventory or assembly changed during host check")
    result.update(status="assembly_glibc_floor_preflight_passed", inventory=pinned,
                  assembly=validated, loaders=loaders, source=record(__file__))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--assembly-manifest-sha256", required=True)
    parser.add_argument("--abi-inventory", type=Path, required=True)
    parser.add_argument("--abi-inventory-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = check_host(args.root, args.abi_inventory, args.abi_inventory_sha256, args.assembly_manifest_sha256)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
