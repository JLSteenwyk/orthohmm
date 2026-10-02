"""Package committed scientific/workflow source separately; verify offline."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess

SCIENTIFIC_REVISION = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
RUNNER = "benchmark_tools/bundle_publication_source.py"
GUIDE = "benchmark_tools/PUBLICATION_SOURCE_COMPONENT.md"
SUPPORT_PINS = {
    "benchmark_tools/results/orthobench_factorial_prepared_20260916.json":
        "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382",
    "benchmark_tools/results/orthobench_paired_uncertainty_20260916.json":
        "660ead29c5b6ac0b8278cd1e62cdcdb0a513db81dda317e634b805e661d70ba9",
    "benchmark_tools/results/integrated_orthobench_data_20260927.json":
        "fcae062525eec61a11de876ae798acffc0f2fe9614466c8ce339a6c214666061",
}
PROFILES = {"source-only", "orthobench-inputs"}


def identity(content):
    return dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def relative(name):
    if not isinstance(name, str):
        raise ValueError("Require a source path string")
    path = PurePosixPath(name)
    if (not path.parts or path.is_absolute()
            or ".." in path.parts or str(path) != name):
        raise ValueError("Unsafe or noncanonical source path")
    return name


def selected(name, component, profile="source-only"):
    if profile not in PROFILES:
        raise ValueError("Unknown source profile")
    if component == "scientific":
        return name in {"LICENSE.md", "README.md", "requirements.txt", "setup.py"} or name.startswith("orthohmm/")
    if component == "workflow":
        return (name in {"LICENSE.md", GUIDE, "benchmark_tools/PUBLICATION_REPRODUCTION.md"}
                or re.fullmatch(r"benchmark_tools/[^/]+\.py", name) is not None
                or re.fullmatch(r"tests/unit/[^/]+\.py", name) is not None
                or profile == "orthobench-inputs" and name in SUPPORT_PINS)
    raise ValueError("Unknown source component")


def required(component, profile):
    if component == "scientific":
        return {"LICENSE.md", "setup.py", "orthohmm/version.py"}
    names = {"LICENSE.md", RUNNER, GUIDE}
    if profile == "orthobench-inputs":
        names |= {*SUPPORT_PINS, "benchmark_tools/verify_orthobench_acquisition.py",
                  "benchmark_tools/rebind_orthobench_data.py"}
    return names


def inventory(repo, revision, component, profile="source-only"):
    data = subprocess.check_output(["git", "-C", str(repo), "ls-tree", "-rz", revision])
    result = {}
    for entry in data.split(b"\0"):
        if not entry:
            continue
        header, raw_name = entry.split(b"\t", 1)
        name = relative(raw_name.decode())
        if not selected(name, component, profile):
            continue
        mode, kind, blob = header.decode().split()
        if kind != "blob" or mode not in {"100644", "100755"} or name in result:
            raise ValueError("Require distinct regular committed source blobs")
        content = subprocess.check_output(["git", "-C", str(repo), "cat-file", "blob", blob])
        if name in SUPPORT_PINS and identity(content)["sha256"] != SUPPORT_PINS[name]:
            raise ValueError("Frozen OrthoBench support manifest differs")
        result[name] = (content, int(mode, 8) & 0o777, blob)
    if not required(component, profile) <= result.keys():
        raise ValueError("Source component is incomplete")
    return result


def verify(directory, manifest_sha):
    directory = Path(directory).resolve(strict=True)
    manifest_path = directory / "SOURCE_INDEX.json"
    if manifest_path.is_symlink():
        raise ValueError("Source index must not be a symlink")
    content = manifest_path.read_bytes()
    if identity(content)["sha256"] != manifest_sha:
        raise ValueError("Source index digest differs from external anchor")
    manifest = json.loads(content)
    if manifest["schema"] == "publication_source_components_v1":
        if "profile" in manifest:
            raise ValueError("Historical source schema cannot override profile")
        profile = "source-only"
    elif manifest["schema"] == "publication_source_components_v2":
        profile = manifest.get("profile")
        if profile != "orthobench-inputs":
            raise ValueError("Unsupported acquisition-support profile")
    else:
        raise ValueError("Unknown source schema")
    if (
            manifest["scientific_revision"] != SCIENTIFIC_REVISION
            or not re.fullmatch(r"[0-9a-f]{40}", manifest["workflow_revision"])
            or manifest["publication_ready"] is not False
            or manifest["redistribution_clearance"] is not False):
        raise ValueError("Source component scope differs")
    seen, counts, syntax, total = set(), {"scientific": 0, "workflow": 0}, 0, 0
    for row in manifest["files"]:
        name = relative(row["path"])
        component, source = name.split("/", 1)
        expected_revision = SCIENTIFIC_REVISION if component == "scientific" else manifest["workflow_revision"]
        if (name in seen or not selected(source, component, profile) or row["git_path"] != source
                or row["git_revision"] != expected_revision or row["mode"] not in (0o644, 0o755)
                or not re.fullmatch(r"[0-9a-f]{40}", row["git_blob"])):
            raise ValueError("Source index contains invalid mapping or duplicate")
        path = directory / name
        if (path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(directory)
                or path.stat().st_mode & 0o777 != row["mode"]):
            raise ValueError("Source payload type, location or mode differs")
        payload = path.read_bytes()
        if identity(payload) != {key: row[key] for key in ("bytes", "sha256")}:
            raise ValueError("Source payload identity differs")
        if source in SUPPORT_PINS and identity(payload)["sha256"] != SUPPORT_PINS[source]:
            raise ValueError("Frozen OrthoBench support manifest differs")
        if source.endswith(".py"):
            compile(payload, name, "exec")
            syntax += 1
        counts[component] += 1
        total += len(payload)
        seen.add(name)
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
    if actual != seen | {"SOURCE_INDEX.json"}:
        raise ValueError("Extra or missing source payloads")
    required_paths = {component + "/" + name for component in counts
                      for name in required(component, profile)}
    if not required_paths <= seen or not all(counts.values()):
        raise ValueError("Missing source component support")
    return dict(status="publication_source_components_verified", files=len(seen), components=counts,
                payload_bytes=total, python_files_syntax_checked=syntax,
                manifest=identity(content), scientific_revision=SCIENTIFIC_REVISION,
                workflow_revision=manifest["workflow_revision"], executable_benchmark_reproduced=False,
                redistribution_clearance=False, publication_ready=False)


def build(repo, revision, output, profile="source-only"):
    if profile not in PROFILES:
        raise ValueError("Unknown source profile")
    repo, output = Path(repo).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", revision + "^{commit}"], text=True).strip()
    payloads, rows = {}, []
    for component, source_revision in (("scientific", SCIENTIFIC_REVISION), ("workflow", commit)):
        for source, (content, mode, blob) in sorted(inventory(repo, source_revision, component, profile).items()):
            name = component + "/" + source
            payloads[name] = (content, mode)
            rows.append(dict(path=name, git_path=source, git_revision=source_revision,
                             git_blob=blob, mode=mode, **identity(content)))
    manifest = dict(schema="publication_source_components_v1", scientific_revision=SCIENTIFIC_REVISION,
                    workflow_revision=commit, files=rows, publication_ready=False,
                    redistribution_clearance=False,
                    exclusions=["Datasets and reference files", "Native predictions and scoring outputs",
                        "Third-party wheels, binaries and source archives", "Runtime/OS images",
                        "Result receipts, plans, figures and manuscript assets", "Test sample datasets"],
                    limitations=["Distinct scientific/workflow revisions are intentional; neither is silently upgraded.",
                        "Integrity and syntax checks do not execute imports, install dependencies or reproduce inference.",
                        "Historical absolute paths, remote hosts and unavailable assets in source are not rewritten or authorized.",
                        "Project license copies do not establish third-party attribution or complete release clearance."])
    if profile == "orthobench-inputs":
        manifest.update(schema="publication_source_components_v2", profile=profile)
        manifest["exclusions"][4] = "Result receipts/plans other than three fixed acquisition-support manifests; figures and manuscript assets"
        manifest["limitations"].append(
            "Support manifests preserve historical provenance paths; acquisition/rebinding must use separately supplied local inputs. No raw data or native runtime is included.")
    output.mkdir(parents=True, exist_ok=False)
    for name, (content, mode) in payloads.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(content)
        path.chmod(mode)
    with (output / "SOURCE_INDEX.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    digest = identity((output / "SOURCE_INDEX.json").read_bytes())["sha256"]
    return verify(output, digest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--repo", type=Path, required=True)
    builder.add_argument("--revision", required=True)
    builder.add_argument("--output", type=Path, required=True)
    builder.add_argument("--profile", choices=sorted(PROFILES), default="source-only")
    verifier = commands.add_parser("verify")
    verifier.add_argument("directory", type=Path)
    verifier.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    result = build(args.repo, args.revision, args.output, args.profile) if args.command == "build" else verify(args.directory, args.manifest_sha256)
    print(json.dumps(result, indent=2, sort_keys=True))
