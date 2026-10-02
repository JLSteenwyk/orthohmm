"""Assemble committed manuscript, source and arithmetic into a local handoff candidate."""

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import subprocess

RUNNER = "benchmark_tools/bundle_publication_handoff.py"
GUIDE = "benchmark_tools/PUBLICATION_HANDOFF.md"
REVIEW_REVISION = "84d4e9bdfc026c10b6cdb170d4467c85daaa3e34"
LEDGER_REVISION = "140b6ef3f78f5f7e3b11396a49c571ed8648a387"
STAGES = dict(render="benchmark_tools/results/publication_main_render_20261002_v3.json",
    print="benchmark_tools/results/publication_main_print_20261002_v3/print.json",
    review="benchmark_tools/results/publication_main_pdf_review_20261002_v3/report.json")
EXTRAS = {
    "arithmetic/reproduce_ygob_validation.py": ("benchmark_tools/reproduce_ygob_validation.py", "d20ab936ed74f031f61dedfaf8cdd69135ef401f1186116d0a97b5ed90a33d16"),
    "arithmetic/reproduce_simulation_panels.py": ("benchmark_tools/reproduce_simulation_panels.py", "f475cb322346a23ae2afeb19a898a9a79fd786645e5db24ce7b3d07221b014b6"),
    "arithmetic/requirements-ygob-arithmetic.txt": ("benchmark_tools/results/ygob_replay_requirements_20261002.txt", "a23aaf94250a9f53031a592980142245a848e2372b6c6c3093e8260b129265b8"),
    "arithmetic/ygob_frozen_results_20260916.json": ("benchmark_tools/results/ygob_frozen_results_20260916.json", "3927d81b1fd851eef3c8ce1343679558ef4f93bcfe2e643435e369587a44a4ef"),
    "arithmetic/ygob_sufficient_counts_20261002.json.gz": ("benchmark_tools/results/ygob_sufficient_counts_20261002.json.gz", "62885bb85eb357d97e55d053134d262e1ebb1e2b088a37df96986f98f8620d66"),
    "arithmetic/simulation_fixed_native_results_20260916.json": ("benchmark_tools/results/simulation_fixed_native_results_20260916.json", "305dcb1dde0c0f57d8b390b0f00cc6148d95103a7efb98bd9f6dce37b96e64a8"),
    "arithmetic/simulation_variable_native_results_20260916.json": ("benchmark_tools/results/simulation_variable_native_results_20260916.json", "cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f"),
    "arithmetic/YGOB_ARITHMETIC_REPLAY.md": ("benchmark_tools/YGOB_ARITHMETIC_REPLAY.md", None),
    "arithmetic/SIMULATION_ARITHMETIC_REPLAY.md": ("benchmark_tools/SIMULATION_ARITHMETIC_REPLAY.md", None),
    "comparison/scores.tsv": ("benchmark_tools/results/current_benchmark_scores_20260926_v2/scores.tsv", "98a38c0e819eefb23258fa85139fdc7a92a43b0e070c8cf8e7035f93381d01ca"),
    "comparison/scores.md": ("benchmark_tools/results/current_benchmark_scores_20260926_v2/scores.md", "9f21cfcb2ea3d237a19faa5a431b3db6d4c45d3358ca9cb22d123b08ff284696"),
    "comparison/provenance.json": ("benchmark_tools/results/current_benchmark_scores_20260926_v2/manifest.json", "8d56dae74721c4530ab2b86c5008b2d0ebac07878e574fd85b1ad84e4b908f07"),
    "README.md": (GUIDE, None),
    "bundle_publication_handoff.py": (RUNNER, None),
    "LICENSE.md": ("LICENSE.md", None),
}
COMPONENTS = dict(source=("SOURCE_INDEX.json", "workflow/benchmark_tools/bundle_publication_source.py"),
                  manuscript=("REVIEW_INDEX.json", "benchmark_tools/bundle_publication_review.py"))
NATIVE_EXTRAS = {
    "runtime/PUBLICATION_RUNTIME_ASSEMBLY_20261002.md": ("benchmark_tools/results/PUBLICATION_RUNTIME_ASSEMBLY_20261002.md", None),
    "runtime/publication_runtime_assembly_20261002.json": ("benchmark_tools/results/publication_runtime_assembly_20261002.json", "aa579ae061cf9a7a7f373902b4b1342495cefe2c5ecc1f10eb894fc7f2d45991"),
    "runtime/publication_runtime_assembly_validation_20261002.json": ("benchmark_tools/results/publication_runtime_assembly_validation_20261002.json", "df975e27c3793cf49578d76dca38b975b98022e4d244fffa590a09f5352ec2cc"),
}
PROFILES = {"native-preparation", "native-build"}


def selection(profile):
    if profile not in PROFILES:
        raise ValueError("Unknown handoff source profile")
    return EXTRAS if profile == "native-preparation" else {**EXTRAS, **NATIVE_EXTRAS}


def identity(content):
    return dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def relative(name):
    if not isinstance(name, str):
        raise ValueError("Require a relative path string")
    path = PurePosixPath(name)
    if not path.parts or path.is_absolute() or ".." in path.parts or str(path) != name or "\\" in name:
        raise ValueError("Unsafe handoff path")
    return name


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    # Do not add bytecode caches to the immutable, already verified component.
    exec(compile(Path(path).read_bytes(), str(path), "exec"), module.__dict__)
    return module


def check_table(content):
    reader = csv.DictReader(io.StringIO(content.decode()), delimiter="\t")
    metrics = ["OrthoBench", "GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS",
               "QfO_secondary_mean", "ThreeKingdoms"]
    if reader.fieldnames != ["Method", *metrics]:
        raise ValueError("Unexpected comparison table columns")
    rows = list(reader)
    if len(rows) != 8 or len({r["Method"] for r in rows}) != 8:
        raise ValueError("Require eight distinct comparator rows")
    for row in rows:
        values = {key: float(row[key]) for key in metrics}
        if not all(math.isfinite(v) and 0 <= v <= 1 for v in values.values()):
            raise ValueError("Nonfinite/out-of-range comparison value")
        mean = sum(values[k] for k in ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")) / 6
        if abs(mean - values["QfO_secondary_mean"]) > 1e-12:
            raise ValueError("Project-defined secondary mean differs")
    return dict(methods=8, numeric_cells=72, raw_scoring_repeated=False, secondary_mean_official=False)


def check_runtime_evidence(directory):
    execution_path = directory / "runtime/publication_runtime_assembly_20261002.json"
    validation_path = directory / "runtime/publication_runtime_assembly_validation_20261002.json"
    execution = json.loads(execution_path.read_bytes())
    validation = json.loads(validation_path.read_bytes())
    if (execution["status"] != "assembled_runtime_integrated_fixture_verified"
            or execution["scientific_revision"] != "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
            or execution["workflow_revision"] != "7bedb195536faf3dcd3bc6079c82252706cde3ed"
            or any(execution[key] is not False for key in ("controlled_timing", "publication_ready",
                "redistribution_clearance", "security_clearance", "new_biological_validation", "full_orthobench_rerun",
                "native_execution_repeated", "retry"))
            or execution["base_unchanged"] is not True or execution["base_site_files"] != 883
            or execution["installed_matched_files"] != 5754
            or execution["installed_wheel_counts"] != {"inference": 11, "reader": 5}
            or execution["original_failure_preserved"] is not True
            or len(execution["executor_outcomes"]) != 10
            or any(row["returncode"] != 0 for row in execution["executor_outcomes"])
            or validation["status"] != "assembled_runtime_receipt_independently_verified"
            or {key: validation["receipt"][key] for key in ("bytes", "sha256")} != identity(execution_path.read_bytes())
            or validation["checked_file_identities"] != 160
            or validation["source_external_command_anchor"] is not True
            or validation["assembly_external_command_anchor"] is not True):
        raise ValueError("Retained runtime evidence scope differs")
    return dict(recorded_executor_stages=10, recorded_installed_payload_files=5754,
                recorded_base_site_files=883, original_failure_preserved=True,
                raw_path_readback_repeated=False, native_inference_repeated=False)


def verify(directory, manifest_sha):
    directory = Path(directory).resolve(strict=True)
    index_path = directory / "HANDOFF_INDEX.json"
    if index_path.is_symlink() or identity(index_path.read_bytes())["sha256"] != manifest_sha:
        raise ValueError("Handoff index differs from external anchor")
    index = json.loads(index_path.read_bytes())
    if index["schema"] == "publication_handoff_candidate_v1":
        profile = "native-preparation"
        if "source_profile" in index:
            raise ValueError("Legacy handoff cannot select a new source profile")
    elif index["schema"] == "publication_handoff_candidate_v2" and index.get("source_profile") == "native-build":
        profile = "native-build"
    else:
        raise ValueError("Handoff schema/source profile differs")
    extras = selection(profile)
    if (any(index[key] is not False for key in ("publication_ready", "redistribution_clearance", "public_release_uploaded"))
            or not re.fullmatch(r"[0-9a-f]{40}", index["workflow_revision"])
            or set(index["components"]) != set(COMPONENTS)):
        raise ValueError("Handoff scope differs")
    seen, total = set(), 0
    for row in index["files"]:
        name = relative(row["path"])
        path = directory / name
        if (name in seen or path.is_symlink() or not path.is_file()
                or not path.resolve().is_relative_to(directory) or row["mode"] not in (0o644, 0o755)
                or path.stat().st_mode & 0o777 != row["mode"]
                or type(row["bytes"]) is not int or row["bytes"] < 0
                or not isinstance(row["sha256"], str) or not re.fullmatch(r"[0-9a-f]{64}", row["sha256"])):
            raise ValueError("Invalid handoff inventory")
        payload = path.read_bytes()
        if identity(payload) != {k: row[k] for k in ("bytes", "sha256")}:
            raise ValueError("Handoff payload differs")
        total += len(payload)
        seen.add(name)
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
    if actual != seen | {"HANDOFF_INDEX.json"}:
        raise ValueError("Extra/missing handoff payloads")
    if set(index["extra_sources"]) != set(extras):
        raise ValueError("Handoff arithmetic/support selection differs")
    for name, (git_path, digest) in extras.items():
        row = index["extra_sources"][name]
        if (name not in seen or row["git_path"] != git_path or row["git_revision"] != index["workflow_revision"]
                or not re.fullmatch(r"[0-9a-f]{40}", row["git_blob"])
                or digest is not None and identity((directory / name).read_bytes())["sha256"] != digest):
            raise ValueError("Committed handoff support mapping/pin differs")
    # Verify all bytes before loading the two included, externally anchored stdlib verifiers.
    results, component_files = {}, set()
    for role, (index_name, verifier) in COMPONENTS.items():
        component = index["components"][role]
        if component["index"] != role + "/" + index_name:
            raise ValueError("Component index location differs")
        module = load(directory / role / verifier, "handoff_" + role)
        result = module.verify(directory / role, component["sha256"])
        if result != component["verified"] or result["publication_ready"] is not False:
            raise ValueError("Component semantic verification differs")
        results[role] = result
        child = json.loads((directory / role / index_name).read_bytes())
        component_files |= {role + "/" + row["path"] for row in child["files"]} | {role + "/" + index_name}
        if role == "source" and child.get("profile") != profile:
            raise ValueError("Source profile differs from selected handoff format")
    if seen != component_files | set(extras):
        raise ValueError("Unaccounted handoff payloads")
    if results["manuscript"]["page_count"] != 9:
        raise ValueError("Require the retained nine-page manuscript")
    table = check_table((directory / "comparison/scores.tsv").read_bytes())
    result = dict(status="publication_handoff_candidate_verified", files=len(seen), payload_bytes=total,
        manifest=identity(index_path.read_bytes()), components=results, comparison=table,
        publication_ready=False, public_release_uploaded=False, redistribution_clearance=False,
        native_inference_reproduced=False, numerical_replay_executed=False)
    if profile == "native-build":
        result.update(source_profile=profile, runtime_evidence_included=True, runtime_payload_delivered=False,
                      runtime_evidence=check_runtime_evidence(directory))
    return result


def build(repo, revision, output, source_profile="native-preparation"):
    repo, output = Path(repo).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    selected = selection(source_profile)
    from benchmark_tools import bundle_publication_source as source
    from benchmark_tools import bundle_publication_review as review
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", revision + "^{commit}"], text=True).strip()
    extras, mappings = {}, {}
    for name, (git_path, digest) in selected.items():
        content, mode, blob = review.committed(repo, commit, git_path)
        if digest is not None and identity(content)["sha256"] != digest:
            raise ValueError("Frozen handoff input differs")
        extras[name] = (content, mode)
        mappings[name] = dict(git_path=git_path, git_revision=commit, git_blob=blob)
    output.mkdir(parents=True, exist_ok=False)
    try:
        results = dict(source=source.build(repo, commit, output / "source", source_profile),
            manuscript=review.build(repo, REVIEW_REVISION, LEDGER_REVISION, commit,
                                    output / "manuscript", stages=STAGES))
        # Component builders inherit the host umask for their generated indexes.
        for role, (index_name, _) in COMPONENTS.items():
            (output / role / index_name).chmod(0o644)
        for name, (content, mode) in extras.items():
            target = output / name
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as stream:
                stream.write(content)
            target.chmod(mode)
        files = []
        for path in sorted(output.rglob("*")):
            if path.is_symlink():
                raise ValueError("Unexpected exported symlink")
            if path.is_file():
                files.append(dict(path=path.relative_to(output).as_posix(), mode=path.stat().st_mode & 0o777,
                                  **identity(path.read_bytes())))
        index = dict(schema="publication_handoff_candidate_v1", workflow_revision=commit,
            components={role: dict(index=role + "/" + COMPONENTS[role][0],
                sha256=result["manifest"]["sha256"], verified=result) for role, result in results.items()},
            extra_sources=mappings, files=files, publication_ready=False, redistribution_clearance=False,
            public_release_uploaded=False, limitations=[
                "Local working candidate, not a complete study runtime, cleared public distribution, submission or DOI.",
                "Manuscript is the dated retained nine-page snapshot; source and numerical support have separate revision identities.",
                "Main-text direct assets are included, not every linked document's transitive files or external URLs.",
                "Raw datasets, native predictions, dependency wheels/tools, Conda bootstrap and OS runtime are excluded.",
                "Verification does not execute arithmetic, infer orthology, establish independent accuracy or admit timing.",
                "Historical paths in metadata are provenance, not required reads for verification or the two standalone arithmetic replays."])
        if source_profile == "native-build":
            index.update(schema="publication_handoff_candidate_v2", source_profile=source_profile)
            index["limitations"].append(
                "Native-build adds frozen setup overlay, seven acquisition/runtime support documents and pinned integration evidence; raw/native/runtime payloads remain excluded.")
        with (output / "HANDOFF_INDEX.json").open("x") as stream:
            json.dump(index, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        (output / "HANDOFF_INDEX.json").chmod(0o644)
        return verify(output, identity((output / "HANDOFF_INDEX.json").read_bytes())["sha256"])
    except Exception as error:
        with (output / "FAILED_BUILD.json").open("x") as stream:
            json.dump(dict(status="handoff_build_failed", retry=False, error=str(error)), stream)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--repo", type=Path, required=True)
    builder.add_argument("--revision", required=True)
    builder.add_argument("--output", type=Path, required=True)
    builder.add_argument("--source-profile", choices=sorted(PROFILES), default="native-preparation")
    verifier = commands.add_parser("verify")
    verifier.add_argument("directory", type=Path)
    verifier.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    if args.command == "build":
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        result = build(args.repo, args.revision, args.output, args.source_profile)
    else:
        result = verify(args.directory, args.manifest_sha256)
    print(json.dumps(result, indent=2, sort_keys=True))
