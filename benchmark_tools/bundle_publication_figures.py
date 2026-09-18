"""Export committed direct figure evidence, or verify a relocated export."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess

AUDIT = "benchmark_tools/results/publication_figure_integrity_20260918_v2.json"
RUNNER = "benchmark_tools/bundle_publication_figures.py"
HELPER = "benchmarks/work/publication_method_native_v2/benchmark_tools/replay_high_sensitivity.py"
HELPER_COMMIT = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
HELPER_SOURCE = "benchmark_tools/replay_high_sensitivity.py"
HELPER_SHA = "852c3e4fc1a53de7e6046aa78324da587376c0df6a0db1cd8265af65c75bea0f"


def identity(content):
    return {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def safe_relative(name):
    path = PurePosixPath(name)
    if path.is_absolute() or not path.parts or ".." in path.parts or str(path) != name:
        raise ValueError("Unsafe or noncanonical relative path")
    return name


def references(value):
    if isinstance(value, dict):
        if "sha256" in value and "path" in value:
            if set(value) != {"path", "bytes", "sha256"}:
                raise ValueError("Unsupported historical file record")
            yield value
        else:
            for child in value.values():
                yield from references(child)
    elif isinstance(value, list):
        for child in value:
            yield from references(child)


def git_bytes(repo, commit, name):
    safe_relative(name)
    return subprocess.check_output(["git", "show", f"{commit}:{name}"], cwd=repo)


def verify(directory):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / "bundle.json").read_text())
    if (manifest["schema_version"] != 1 or manifest["scope"] != "direct_figure_evidence_only"
            or manifest["publication_ready"] is not False):
        raise ValueError("Unsupported bundle scope")
    files = {}
    for item in manifest["files"]:
        name = safe_relative(item["path"])
        if name in files:
            raise ValueError("Duplicate bundle path")
        path = directory / name
        if path.is_symlink() or not path.resolve().is_relative_to(directory):
            raise ValueError("Symlink or escaping bundle path")
        if identity(path.read_bytes()) != {k: item[k] for k in ("bytes", "sha256")}:
            raise ValueError("Bundle file identity mismatch")
        files[name] = item
    actual = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
    if actual != set(files) | {"bundle.json"}:
        raise ValueError("Unexpected or missing bundle files")
    used = set()
    panels = set()
    outputs = 0
    for panel in manifest["panels"]:
        if panel["name"] in panels:
            raise ValueError("Duplicate panel")
        panels.add(panel["name"])
        source = panel["manifest"]
        if source not in files:
            raise ValueError("Unlisted figure manifest")
        used.add(source)
        data = json.loads((directory / source).read_text())
        mapping = panel["relocations"]
        expected = {}
        for row in references(data):
            if row["path"] in expected and expected[row["path"]] != row:
                raise ValueError("Conflicting historical records")
            expected[row["path"]] = row
        if set(mapping) != set(expected):
            raise ValueError("Incomplete or extra relocation mapping")
        for historical, relative in mapping.items():
            if relative not in files or any(files[relative][k] != expected[historical][k] for k in ("bytes", "sha256")):
                raise ValueError("Historical evidence identity mismatch")
            used.add(relative)
        recorded_outputs = data["outputs"]
        paths = [r["path"] for r in recorded_outputs]
        if not paths or len(set(paths)) != len(paths) or not {".png", ".pdf", ".svg"} <= {PurePosixPath(p).suffix for p in paths}:
            raise ValueError("Invalid figure output inventory")
        if len(paths) != panel["output_count"]:
            raise ValueError("Output count mismatch")
        outputs += len(paths)
    if used | set(manifest["support_files"]) != set(files):
        raise ValueError("Unaccounted bundle files")
    return {"status": "relocated_direct_figure_evidence_verified", "panels": len(panels),
            "outputs": outputs, "files": len(files), "bytes": sum(r["bytes"] for r in files.values()),
            "bundle_manifest": identity((directory / "bundle.json").read_bytes()),
            "publication_ready": False}


def build(repo, revision, output, audit_path=AUDIT):
    repo, output = Path(repo).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=repo, text=True).strip()
    audit = json.loads(git_bytes(repo, commit, audit_path))
    if audit["status"] != "retained_figure_bytes_verified" or not audit["panels"]:
        raise ValueError("Require successful committed figure audit")
    payloads, entries, panels = {}, {}, []

    def add(name, expected=None):
        safe_relative(name)
        revision, source = (HELPER_COMMIT, HELPER_SOURCE) if name == HELPER else (commit, name)
        if name not in payloads:
            content = git_bytes(repo, revision, source)
            if name == HELPER and identity(content)["sha256"] != HELPER_SHA:
                raise ValueError("Frozen helper mismatch")
            payloads[name] = content
            entries[name] = {"path": name, **identity(content), "git_commit": revision, "git_path": source}
        if expected and any(entries[name][key] != expected[key] for key in ("bytes", "sha256")):
            raise ValueError("Committed bytes differ from retained figure evidence")

    support = [audit_path, "LICENSE.md", RUNNER]
    for name in support:
        add(name)
    for panel in audit["panels"]:
        name = safe_relative(panel["panel"])
        manifest_path = f"benchmark_tools/results/{name}/manifest.json"
        if panel["status"] != "all_recorded_bytes_match":
            raise ValueError("Unverified panel")
        add(manifest_path, panel["manifest"])
        relocations = {}
        for item in panel["files"]:
            if item["status"] != "matches" or item["actual"] != item["expected"]:
                raise ValueError("Unverified historical dependency")
            relative = item["repository_relative_path"]
            if not relative:
                raise ValueError("External dependency has no export mapping")
            add(relative, item["expected"])
            historical = item["expected"]["path"]
            if historical in relocations and relocations[historical] != relative:
                raise ValueError("Ambiguous historical relocation")
            relocations[historical] = relative
        panels.append({"name": name, "manifest": manifest_path, "relocations": relocations,
                       "output_count": panel["output_count"]})
    manifest = {"schema_version": 1, "scope": "direct_figure_evidence_only", "publication_ready": False,
                "source_commit": commit, "files": [entries[k] for k in sorted(entries)],
                "support_files": support, "panels": panels,
                "limitations": ["Historical absolute paths are provenance; use the explicit relocation mappings.",
                    "No native inference, scoring or plotting is executed by this verification.",
                    "Imports, compiled dependencies and transitive data are not included.",
                    "Repository licensing does not establish third-party data redistribution permission.",
                    "Preserves diagnostic and descriptive figures, not new scientific admissions.",
                    "Not the complete publication archive or a release-readiness claim."]}
    output.mkdir(parents=True, exist_ok=False)
    for name, content in payloads.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(content)
    with (output / "bundle.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return verify(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--repo", type=Path, required=True)
    builder.add_argument("--revision", required=True)
    builder.add_argument("--output", type=Path, required=True)
    verifier = commands.add_parser("verify")
    verifier.add_argument("directory", type=Path)
    args = parser.parse_args()
    result = build(args.repo, args.revision, args.output) if args.command == "build" else verify(args.directory)
    print(json.dumps(result, indent=2, sort_keys=True))
