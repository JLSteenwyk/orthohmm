"""Inspect pinned wheel ELF dependencies without loading or executing their code."""

import argparse
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tempfile
import zipfile

from benchmark_tools.export_dependency_notices import safe_member, select_wheels
from benchmark_tools.inventory_dependency_notices import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


def dynamic_tags(text):
    values = {key: [] for key in ("NEEDED", "SONAME", "RPATH", "RUNPATH")}
    for line in text.splitlines():
        for key in values:
            if "(" + key + ")" not in line:
                continue
            match = re.fullmatch(r"\s*0x[0-9a-fA-F]+\s+\(" + key + r"\)\s+[^\[]+\[([^\]]*)\]\s*", line)
            if not match:
                raise ValueError("Unrecognized dynamic string tag: " + line)
            values[key].append(match[1])
    for key in ("SONAME", "RPATH", "RUNPATH"):
        if len(values[key]) > 1:
            raise ValueError("Duplicate " + key)
    return values


def abi_requirements(text):
    fields = ("Class", "Data", "OS/ABI", "ABI Version", "Type", "Machine", "Flags", "Number of program headers")
    header = {}
    for line in text.split("Program Headers:", 1)[0].splitlines():
        key, separator, value = line.strip().partition(":")
        if separator and key in fields:
            if key in header or not value.strip():
                raise ValueError("Duplicate/empty ELF header field")
            header[key] = value.strip()
    if set(header) != set(fields):
        raise ValueError("Incomplete ELF ABI header")
    if not header["Number of program headers"].isdigit():
        raise ValueError("Unsupported program-header count")
    if (int(header["Number of program headers"]) and "Program Headers:" not in text
            or not int(header["Number of program headers"]) and "There are no program headers" not in text):
        raise ValueError("Missing program-header evidence")
    program_rows = re.findall(r"^\s*[A-Z][A-Z0-9_]*\s+0x[0-9a-fA-F]+\s+0x", text, re.M)
    if len(program_rows) != int(header["Number of program headers"]):
        raise ValueError("Incomplete program-header table")
    interpreters = re.findall(r"^\s*\[Requesting program interpreter: ([^\]\n]+)\]\s*$", text, re.M)
    segments = re.findall(r"^\s*INTERP\s+0x", text, re.M)
    if len(interpreters) != len(segments) or len(interpreters) > 1:
        raise ValueError("Incomplete/duplicate program interpreter")
    requirements = []
    active = False
    declared = None
    expected_names = None
    current = None
    for line in text.splitlines():
        if line.startswith("Version needs section "):
            if declared is not None:
                raise ValueError("Duplicate version needs section")
            match = re.fullmatch(r"Version needs section '[^']+' contains (\d+) entr(?:y|ies):", line)
            if not match:
                raise ValueError("Unrecognized version needs section")
            declared = int(match[1])
            active = True
            continue
        if line.startswith("Version "):
            active = False
        if not active or not line.strip():
            continue
        if re.fullmatch(r"\s*Addr: 0x[0-9a-fA-F]+\s+Offset: 0x[0-9a-fA-F]+\s+Link: \d+ \([^)]*\)\s*", line):
            continue
        file_match = re.fullmatch(r"\s*(?:0x)?[0-9a-fA-F]+: Version: 1\s+File: (\S+)\s+Cnt: (\d+)\s*", line)
        if file_match:
            if current is not None and len(current["versions"]) != expected_names:
                raise ValueError("Incomplete version requirements")
            current = dict(library=file_match[1], versions=[])
            expected_names = int(file_match[2])
            requirements.append(current)
            continue
        name_match = re.fullmatch(r"\s*(?:0x)?[0-9a-fA-F]+:\s+Name: (\S+)\s+Flags: (.*?)\s+Version: (\d+)\s*", line)
        if name_match and current is not None:
            current["versions"].append(dict(name=name_match[1], flags=name_match[2], index=int(name_match[3])))
            continue
        raise ValueError("Unrecognized version requirement: " + line)
    if ((current is not None and len(current["versions"]) != expected_names)
            or (declared is not None and len(requirements) != declared)):
        raise ValueError("Incomplete version requirements")
    if len({row["library"] for row in requirements}) != len(requirements):
        raise ValueError("Duplicate required library")
    for row in requirements:
        if len({v["name"] for v in row["versions"]}) != len(row["versions"]):
            raise ValueError("Duplicate required version")
    if declared is None and not ("No version information found in this file." in text
            or re.search(r"^Version (?:symbols|definition) section ", text, re.M)):
        raise ValueError("Missing version-information evidence")
    return dict(header=header, interpreter=interpreters[0] if interpreters else None,
                version_requirements=requirements)


def scan_abi(binary, readelf):
    result = subprocess.run([str(readelf), "--wide", "--file-header", "--program-headers", "--version-info", str(binary)],
                            capture_output=True, text=True, check=True,
                            env={"LC_ALL": "C", "PATH": "/usr/bin:/bin"}, timeout=60)
    if result.stderr:
        raise ValueError("readelf ABI diagnostics: " + result.stderr)
    return dict(requirements=abi_requirements(result.stdout), readelf_stdout=result.stdout)


def scan_wheel(wheel, readelf, *, include_abi=False):
    check(wheel["wheel"])
    objects = []
    with zipfile.ZipFile(wheel["wheel"]["path"]) as archive, tempfile.TemporaryDirectory() as temp:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate wheel members")
        for info in archive.infolist():
            if info.is_dir():
                continue
            safe_member(info.filename)
            if stat.S_IFMT(info.external_attr >> 16) not in (0, stat.S_IFREG):
                raise ValueError("Nonregular wheel member")
            with archive.open(info) as source:
                magic = source.read(4)
                if magic != b"\x7fELF":
                    continue
                # Fixed temporary name: archive paths never control extraction paths.
                binary = Path(temp) / "object"
                with binary.open("wb") as dest:
                    dest.write(magic)
                    shutil.copyfileobj(source, dest)
            identity = record(binary)
            command = [str(readelf), "--wide", "--dynamic", str(binary)]
            result = subprocess.run(command, capture_output=True, text=True, check=True,
                                    env={"LC_ALL": "C", "PATH": "/usr/bin:/bin"}, timeout=60)
            if result.stderr:
                raise ValueError("readelf diagnostics: " + result.stderr)
            obj = dict(member=info.filename, bytes=identity["bytes"], sha256=identity["sha256"],
                       tags=dynamic_tags(result.stdout), readelf_stdout=result.stdout)
            if include_abi:
                obj["abi"] = scan_abi(binary, readelf)
            objects.append(obj)
    check(wheel["wheel"])
    return dict(name=wheel["name"], version=wheel["version"], wheel=wheel["wheel"],
                objects=objects, notice_candidates=wheel["notice_candidates"],
                filename_candidates_not_elf=sorted(set(wheel["native_members"]) - {o["member"] for o in objects}))


def candidate_edges(wheels):
    providers = {}
    for wheel in wheels:
        for obj in wheel["objects"]:
            identity = dict(wheel_sha256=wheel["wheel"]["sha256"], member=obj["member"])
            for name in set([PurePosixPath(obj["member"]).name, *obj["tags"]["SONAME"]]):
                providers.setdefault(name, []).append(identity)
    edges = []
    for wheel in wheels:
        for obj in wheel["objects"]:
            for name in obj["tags"]["NEEDED"]:
                edges.append(dict(wheel_sha256=wheel["wheel"]["sha256"], member=obj["member"], needed=name,
                                  selected_wheel_candidates=providers.get(name, [])))
    return edges


def inspect(paths, readelf):
    binary = record(readelf)
    watched, inventories = [], []
    for path in paths:
        watched.append(record(path))
        inventory = json.loads(path.read_text())
        report = inventory["install_report"]
        if audit(Path(report["path"]), report["sha256"]) != inventory:
            raise ValueError("Notice inventory does not reproduce")
        inventories.append(inventory)
    wheels = [scan_wheel(w, readelf) for w in select_wheels(inventories)]
    edges = candidate_edges(wheels)
    for item in [binary, *watched]:
        check(item)
    return dict(status="selected_wheel_elf_inventory", source=record(__file__), inputs=watched,
                readelf=binary, readelf_version=subprocess.check_output([str(readelf), "--version"], text=True),
                wheels=wheels, dependency_edges=edges,
                needs_without_selected_wheel_candidate=sorted({e["needed"] for e in edges if not e["selected_wheel_candidates"]}),
                publication_ready=False, redistribution_clearance=False,
                limitations=["NEEDED and search paths are declarations, not observed runtime resolution.",
                    "Candidates match basename or SONAME across selected wheels, ignoring loader scope and search order.",
                    "Static libraries, dlopen dependencies, external tools, Python and OS libraries need separate review.",
                    "Wheel notice candidates are not component-level license attribution or compatibility judgments.",
                    "Only ELF content is inspected; no wheel code is loaded or executed."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", action="append", type=Path, required=True)
    parser.add_argument("--readelf", type=Path, default=Path("/usr/bin/readelf"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    save(args.output, inspect(args.inventory, args.readelf.resolve()))
