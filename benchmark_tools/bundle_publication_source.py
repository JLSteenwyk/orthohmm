"""Package committed scientific/workflow source separately; verify offline."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess

SCIENTIFIC_REVISION = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
BUILD_REVISION = "6fd6df19daba83ec6467b917988f99e27a95be14"
BUILD_SETUP_SHA = "88120e9d722557337d466a5026c4b238d5a21c9cc213bf35c27d4b480b347193"
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
BASE_SUPPORT_PINS = {
    "benchmark_tools/results/reconstructed_base_fixture_20260927.json":
        "7af55a78192eea22cf978ff331c1c5556965e17862b565793321b14857479f59",
}
BASE_HELPERS = {"benchmark_tools/" + name + ".py" for name in (
    "install_publication_base", "stage_base_archives", "run_integrated_publication_workflow",
    "audit_recovery_install", "audit_frozen_overlay_install", "audit_leiden_recovery_wheel",
    "prepare_frozen_build_overlay", "verify_frozen_source_archive", "prepare_ob_candidate_neighborhood")}
WHEEL_SUPPORT_PINS = {
    "benchmark_tools/results/integrated_wheel_elf_20260927.json":
        "97eaa8cca56db98b337ad161e619006ef7a91433d2c704e9dcc16ceacc5e7dab",
    "benchmark_tools/results/publication_recovery_requirements_20260926.txt":
        "813751262c9405155cc2f1974a2c32187b352f88307e2332e0cec9df35708ae3",
    "benchmark_tools/results/publication_reader_requirements_20260927_v2.txt":
        "53df1cabd91c2fb18179951ebc95166ae49ba843f4015e0f417468d2a8265063",
}
WHEEL_HELPERS = {"benchmark_tools/acquire_publication_wheels.py", "benchmark_tools/acquire_publication_base.py"}
BUILD_HELPERS = {"benchmark_tools/build_publication_project_wheel.py"}
PROFILES = {"source-only", "orthobench-inputs", "native-preparation", "native-wheels", "native-build"}


def support_pins(profile):
    if profile in {"native-wheels", "native-build"}:
        return {**SUPPORT_PINS, **BASE_SUPPORT_PINS, **WHEEL_SUPPORT_PINS}
    if profile == "native-preparation":
        return {**SUPPORT_PINS, **BASE_SUPPORT_PINS}
    return SUPPORT_PINS if profile == "orthobench-inputs" else {}


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
    if component == "build":
        return profile == "native-build" and name == "setup.py"
    if component == "workflow":
        return (name in {"LICENSE.md", GUIDE, "benchmark_tools/PUBLICATION_REPRODUCTION.md"}
                or re.fullmatch(r"benchmark_tools/[^/]+\.py", name) is not None
                or re.fullmatch(r"tests/unit/[^/]+\.py", name) is not None
                or name in support_pins(profile))
    raise ValueError("Unknown source component")


def required(component, profile):
    if component == "scientific":
        return {"LICENSE.md", "setup.py", "orthohmm/version.py"}
    if component == "build":
        return {"setup.py"} if profile == "native-build" else set()
    names = {"LICENSE.md", RUNNER, GUIDE}
    if profile in {"orthobench-inputs", "native-preparation", "native-wheels", "native-build"}:
        names |= {*SUPPORT_PINS, "benchmark_tools/verify_orthobench_acquisition.py",
                  "benchmark_tools/rebind_orthobench_data.py"}
    if profile in {"native-preparation", "native-wheels", "native-build"}:
        names |= {*BASE_SUPPORT_PINS, *BASE_HELPERS}
    if profile in {"native-wheels", "native-build"}:
        names |= {*WHEEL_SUPPORT_PINS, *WHEEL_HELPERS}
    if profile == "native-build":
        names |= BUILD_HELPERS
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
        if component == "build" and identity(content)["sha256"] != BUILD_SETUP_SHA:
            raise ValueError("Frozen setup-only build overlay differs")
        if name in support_pins(profile) and identity(content)["sha256"] != support_pins(profile)[name]:
            raise ValueError("Frozen acquisition/runtime support manifest differs")
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
        if profile not in PROFILES - {"source-only"}:
            raise ValueError("Unsupported acquisition-support profile")
    else:
        raise ValueError("Unknown source schema")
    if profile == "native-build" and manifest.get("build_revision") != BUILD_REVISION:
        raise ValueError("Build overlay revision differs")
    if (
            manifest["scientific_revision"] != SCIENTIFIC_REVISION
            or not re.fullmatch(r"[0-9a-f]{40}", manifest["workflow_revision"])
            or manifest["publication_ready"] is not False
            or manifest["redistribution_clearance"] is not False):
        raise ValueError("Source component scope differs")
    seen, counts, syntax, total = set(), {"scientific": 0, "workflow": 0}, 0, 0
    if profile == "native-build":
        counts["build"] = 0
    for row in manifest["files"]:
        name = relative(row["path"])
        component, source = name.split("/", 1)
        expected_revision = (SCIENTIFIC_REVISION if component == "scientific" else
                             BUILD_REVISION if component == "build" else manifest["workflow_revision"])
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
        if component == "build" and identity(payload)["sha256"] != BUILD_SETUP_SHA:
            raise ValueError("Frozen setup-only build overlay differs")
        if source in support_pins(profile) and identity(payload)["sha256"] != support_pins(profile)[source]:
            raise ValueError("Frozen acquisition/runtime support manifest differs")
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
    result = dict(status="publication_source_components_verified", files=len(seen), components=counts,
                payload_bytes=total, python_files_syntax_checked=syntax,
                manifest=identity(content), scientific_revision=SCIENTIFIC_REVISION,
                workflow_revision=manifest["workflow_revision"], executable_benchmark_reproduced=False,
                redistribution_clearance=False, publication_ready=False)
    if profile == "native-build":
        result["build_revision"] = BUILD_REVISION
    return result


def build(repo, revision, output, profile="source-only"):
    if profile not in PROFILES:
        raise ValueError("Unknown source profile")
    repo, output = Path(repo).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", revision + "^{commit}"], text=True).strip()
    payloads, rows = {}, []
    components = [("scientific", SCIENTIFIC_REVISION), ("workflow", commit)]
    if profile == "native-build":
        components.append(("build", BUILD_REVISION))
    for component, source_revision in components:
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
    if profile in PROFILES - {"source-only"}:
        manifest.update(schema="publication_source_components_v2", profile=profile)
        manifest["exclusions"][4] = "Result receipts/plans other than three fixed acquisition-support manifests; figures and manuscript assets"
        manifest["limitations"].append(
            "Support manifests preserve historical provenance paths; acquisition/rebinding must use separately supplied local inputs. No raw data or native runtime is included.")
        if profile in {"native-preparation", "native-wheels", "native-build"}:
            manifest["exclusions"][4] = "Result receipts/plans other than four fixed acquisition/runtime support documents; figures and manuscript assets"
            manifest["limitations"].append(
                "Native-preparation additionally includes the fixed historical base reconstruction receipt and required controller helpers, not archives, wheels or Conda bootstrap.")
        if profile in {"native-wheels", "native-build"}:
            manifest["exclusions"][4] = "Result receipts/plans other than seven fixed acquisition/runtime/wheel-support documents; figures and manuscript assets"
            manifest["limitations"].append(
                "Native-wheels includes the exact admitted wheel inventory and both historical hash locks; the unpublished project wheel and pip still require separately supplied exact artifacts.")
        if profile == "native-build":
            manifest["build_revision"] = BUILD_REVISION
            manifest["limitations"].append(
                "Native-build additionally carries the frozen setup-only overlay separately; a rebuilt candidate wheel is not silently substituted into historical hash locks or admitted scientific executors.")
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
