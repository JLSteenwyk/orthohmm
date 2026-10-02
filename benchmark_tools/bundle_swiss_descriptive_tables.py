"""Archive retained SwissTrees descriptive arithmetic, not raw-source admission."""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import runpy
import subprocess
import tarfile
import tempfile

RUNNER = "benchmark_tools/bundle_swiss_descriptive_tables.py"
CHECKER = "benchmark_tools/check_swiss_descriptive_tables.py"
CHECKER_SHA = "43f61640b1e8482d2501731a6fcd547a0287878d297cf5a30acd6249f5b0469c"
SCOPE = "retained_swiss_descriptive_table_arithmetic_only"
MAX_BYTES = 4 * 1024 ** 2
FEATURES = {
    "descriptive": "corrected_swiss_sequence_strata_20260918.json",
    "identity": "corrected_swiss_identity_admission_22102.json",
    "fragment": "swiss_historical_fragment_admission_22117.json",
    "duplication": "swiss_duplication_features_v2_20260923.json",
}
README = """# SwissTrees Descriptive Table Component

Eight methods, eighteen reference families and all four retained descriptive
tables. Only Python's standard library is needed. Historical absolute paths
inside reports are provenance, not accessible raw-input dependencies.

After safe transfer, retain the bundle.json digest independently and run:

```sh
python3 -I -B benchmark_tools/bundle_swiss_descriptive_tables.py verify --directory /absolute/component --manifest-sha256 REVIEWED_SHA256
python3 -I -B benchmark_tools/bundle_swiss_descriptive_tables.py reproduce --directory /absolute/component --manifest-sha256 REVIEWED_SHA256 --output /absolute/new-report.json
```

Reproduction checks 208 JSON/TSV rows and 984 logical score/difference cells
using exact rational family precision/recall and their harmonic macro F1.
Missing cells stay missing. Markdown bytes are checked, not regenerated.
No package installation, Git, bootstrap, plotting, native inference, raw-QfO
scoring or original annotation/alignment/tree admission occurs on replay.

The project license does not grant rights to every upstream dataset. This is
not independent biological confirmation, controlled timing, rights clearance,
the complete executable study release or public archival deposition.
"""


def identity(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def data_paths():
    names = ["qfo_recovered_swiss_uncertainty_22178.json"]
    for kind, feature in FEATURES.items():
        names.append(feature)
        names.extend(f"swiss_{kind}_strata_20260926/{name}" for name in ("manifest.json", "scores.tsv", "scores.md"))
    return names


def expected_files():
    return {RUNNER, CHECKER, "LICENSE.md", "README.md", *["results/" + name for name in data_paths()]}


def safe_name(name):
    path = PurePosixPath(name)
    if (not name or path.is_absolute() or ".." in path.parts or str(path) != name
            or "\\" in name or ":" in name):
        raise ValueError("Unsafe component member path")
    return name


def direct(path):
    path = Path(path)
    if not path.is_absolute() or path.resolve() != path or path.is_symlink():
        raise ValueError("Require direct absolute component path")
    return path


def verify(directory, manifest_sha):
    directory = direct(directory)
    raw = direct(directory / "bundle.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != manifest_sha:
        raise ValueError("Reviewed component manifest changed")
    manifest = json.loads(raw)
    if (manifest["schema_version"] != 1 or manifest["scope"] != SCOPE
            or manifest["publication_ready"] is not False or manifest["raw_source_admission"] is not False
            or not re.fullmatch("[0-9a-f]{40}", manifest["source_commit"])):
        raise ValueError("Wrong component scope")
    expected, payloads = expected_files(), {}
    for ref in manifest["files"]:
        name = safe_name(ref["path"])
        if name in payloads or name not in expected:
            raise ValueError("Duplicate or unexpected component member")
        data = direct(directory / name).read_bytes()
        if identity(data) != {k: ref[k] for k in ("bytes", "sha256")}:
            raise ValueError("Changed component bytes")
        payloads[name] = data
    if set(payloads) != expected:
        raise ValueError("Missing component member")
    observed = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
    if observed != expected | {"bundle.json"}:
        raise ValueError("Unexpected component filesystem member")
    if (hashlib.sha256(payloads[CHECKER]).hexdigest() != CHECKER_SHA
            or payloads[RUNNER] != Path(__file__).read_bytes() or payloads["README.md"] != README.encode()):
        raise ValueError("Changed arithmetic checker, runner or guide")
    return dict(status="swiss_descriptive_component_verified", manifest=identity(raw),
        files=len(payloads), bytes=sum(len(v) for v in payloads.values()),
        source_commit=manifest["source_commit"], numerical_reproduction_executed=False,
        raw_source_admission=False, publication_ready=False), payloads


def reproduce(directory, manifest_sha, output):
    directory, output = direct(directory), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.is_relative_to(directory):
        raise ValueError("Reproduction output must not mutate the component")
    checked, _ = verify(directory, manifest_sha)
    namespace = runpy.run_path(str(directory / CHECKER))
    result = namespace["verify"](directory / "results")
    if (result["rows"] != 208 or result["score_and_difference_cells"] != 984
            or result["raw_inputs_revalidated"] is not False or result["bootstrap_intervals_recomputed"] is not False
            or verify(directory, manifest_sha)[0] != checked):
        raise ValueError("Incomplete arithmetic or changed component")
    result.update(component=checked, source_commit=checked["source_commit"],
                  component_numerical_reproduction_executed=True)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


def build(repo, revision, output):
    repo, output = direct(repo), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", revision + "^{commit}"], text=True).strip()
    payloads, entries = {}, []
    sources = {RUNNER: RUNNER, CHECKER: CHECKER, "LICENSE.md": "LICENSE.md",
               **{"results/" + name: "benchmark_tools/results/" + name for name in data_paths()}}
    for target, source in sorted(sources.items()):
        data = subprocess.check_output(["git", "-C", str(repo), "show", commit + ":" + source])
        payloads[target] = data
        entries.append(dict(path=target, **identity(data), git_path=source))
    if payloads[RUNNER] != Path(__file__).read_bytes() or hashlib.sha256(payloads[CHECKER]).hexdigest() != CHECKER_SHA:
        raise ValueError("Executing builder or frozen checker differs from committed source")
    payloads["README.md"] = README.encode()
    entries.append(dict(path="README.md", **identity(payloads["README.md"]), generated_by=RUNNER))
    manifest = dict(schema_version=1, scope=SCOPE, source_commit=commit, raw_source_admission=False,
                    publication_ready=False, files=sorted(entries, key=lambda r: r["path"]))
    output.mkdir(parents=True, exist_ok=False)
    for name, data in sorted(payloads.items()):
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(data)
    with (output / "bundle.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    digest = identity((output / "bundle.json").read_bytes())["sha256"]
    return verify(output, digest)[0]


def archive(directory, manifest_sha, output):
    directory, output = direct(directory), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.is_relative_to(directory):
        raise ValueError("Archive must not mutate component")
    checked, payloads = verify(directory, manifest_sha)
    payloads["bundle.json"] = (directory / "bundle.json").read_bytes()
    with output.open("xb") as stream, gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode="w") as handle:
            for name, data in sorted(payloads.items()):
                info = tarfile.TarInfo(name)
                info.size, info.mode, info.mtime = len(data), 0o644, 0
                handle.addfile(info, io.BytesIO(data))
    if verify(directory, manifest_sha)[0] != checked:
        raise ValueError("Component changed during archiving")
    return dict(status="swiss_descriptive_component_archived", archive=identity(output.read_bytes()),
                members=len(payloads), component=checked, publication_ready=False)


def restore_and_reproduce(archive_path, manifest_sha, output):
    archive_path, output = direct(archive_path), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if archive_path.stat().st_size > MAX_BYTES:
        raise ValueError("Archive exceeds size budget")
    expected = expected_files() | {"bundle.json"}
    seen, total = set(), 0
    with tempfile.TemporaryDirectory(prefix="swiss-tables-restored-") as temporary:
        root = Path(temporary).resolve()
        with tarfile.open(archive_path, "r|gz") as handle:
            for member in handle:
                name = safe_name(member.name)
                total += member.size
                if (name in seen or name not in expected or not member.isfile() or member.mode != 0o644
                        or member.size < 0 or total > MAX_BYTES):
                    raise ValueError("Unexpected, duplicate or oversized archive member")
                seen.add(name)
                data = handle.extractfile(member).read(member.size + 1)
                if len(data) != member.size:
                    raise ValueError("Archive member size differs")
                target = root / name
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open("xb") as stream:
                    stream.write(data)
        if seen != expected:
            raise ValueError("Incomplete archive inventory")
        result = reproduce(root, manifest_sha, output)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--repo", type=Path, required=True)
    builder.add_argument("--revision", required=True)
    builder.add_argument("--output", type=Path, required=True)
    for name in ("verify", "reproduce", "archive", "restore"):
        command = commands.add_parser(name)
        command.add_argument("--manifest-sha256", required=True)
        command.add_argument("--archive" if name == "restore" else "--directory", type=Path, required=True)
        if name != "verify":
            command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "build":
        result = build(args.repo, args.revision, args.output)
    elif args.command == "verify":
        result = verify(args.directory, args.manifest_sha256)[0]
    elif args.command == "restore":
        result = restore_and_reproduce(args.archive, args.manifest_sha256, args.output)
    else:
        result = globals()[args.command](args.directory, args.manifest_sha256, args.output)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
