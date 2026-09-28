"""Read Git objects to diagnose source-version labels; never build native code."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record

ORIGIN = "https://github.com/vtraag/libleidenalg.git"
COMMIT = "7ce12032cd61a38caebc52881ce1f1a1be4f0d28"
FILES = ("CMakeLists.txt", "src/CMakeLists.txt", "etc/cmake/version.cmake",
         "etc/cmake/GetGitRevisionDescription.cmake", ".gitattributes", "vcpkg.json", "LICENSE")


def git(repo, *args):
    command = ["git", "--no-replace-objects", "-C", str(repo), *args]
    done = subprocess.run(command, capture_output=True, timeout=60)
    if done.returncode:
        raise ValueError(f"Git read failed: {args[0]}: {done.stderr.decode(errors='replace')}")
    return done.stdout


def inspect(repo, commit=COMMIT):
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("Require explicit full commit identity")
    if git(repo, "config", "--get", "remote.origin.url").decode().strip() != ORIGIN:
        raise ValueError("Unexpected source repository origin")
    if git(repo, "rev-parse", "--is-shallow-repository").strip() != b"false":
        raise ValueError("Shallow history cannot reproduce the version lookup")
    if git(repo, "rev-parse", "refs/tags/0.12.0^{commit}").decode().strip() != commit:
        raise ValueError("Release tag points to a different source commit")
    git(repo, "fsck", "--full", "--no-reflogs")
    tree = git(repo, "ls-tree", "-r", "-z", "--full-tree", commit)
    entries, names = [], set()
    selected = {}
    for line in tree.split(b"\0"):
        if not line:
            continue
        metadata, name = line.split(b"\t", 1)
        mode, kind, blob = metadata.decode().split()
        name = name.decode()
        if name in names or kind != "blob" or mode not in {"100644", "100755"}:
            raise ValueError("Duplicate path, symlink or submodule in source inventory")
        names.add(name)
        content = git(repo, "cat-file", "blob", blob)
        actual_blob = hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()
        if actual_blob != blob:
            raise ValueError("Source blob identity differs")
        entries.append(dict(path=name, git_blob=blob, mode=mode, bytes=len(content),
                            sha256=hashlib.sha256(content).hexdigest()))
        if name in FILES:
            selected[name] = content.decode()
    if not set(FILES) <= names:
        raise ValueError("Missing required version-build or license sources")
    if names & {"VERSION", "NEXT_VERSION"}:
        raise ValueError("Version override file would change the inspected build branch")
    if ("git_describe(PACKAGE_VERSION)" not in selected["etc/cmake/version.cmake"]
            or 'string(REGEX MATCH "^[^-]+" PACKAGE_VERSION_BASE "${PACKAGE_VERSION}")' not in selected["etc/cmake/version.cmake"]):
        raise ValueError("Inspected version derivation differs; review source before interpretation")
    description = git(repo, "describe", commit).decode().strip()
    with_tags = git(repo, "describe", "--tags", "--exact-match", commit).decode().strip()
    tags = git(repo, "for-each-ref", "--format=%(refname) %(objecttype) %(objectname) %(*objectname)", "refs/tags").decode()
    # This matches the inspected source's prefix operation, not a CMake execution.
    base = description.split("-", 1)[0]
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", base) or with_tags != "0.12.0":
        raise ValueError("Unexpected observed version labels")
    if git(repo, "rev-parse", "refs/tags/0.12.0^{commit}").decode().strip() != commit:
        raise ValueError("Source tag changed during inspection")
    if git(repo, "for-each-ref", "--format=%(refname) %(objecttype) %(objectname) %(*objectname)", "refs/tags").decode() != tags:
        raise ValueError("Source tag inventory changed during inspection")
    return dict(status="source_git_version_mechanism_observed", source=record(__file__),
        inspected_utc=datetime.now(timezone.utc).isoformat(), repository=str(repo), origin=ORIGIN,
        commit=commit, tree=git(repo, "rev-parse", commit+"^{tree}").decode().strip(),
        git_version=git(repo, "--version").decode().strip(), tag_inventory=tags,
        git_describe=description, git_describe_with_tags=with_tags, inferred_package_version_base=base,
        files=entries, selected_sources=selected, wheel_source_identity_established=False,
        cmake_executed=False, native_code_built=False, publication_ready=False,
        limitations=["Read-only Git/source analysis, not a historical wheel-build attestation.",
                     "The derived version prefix is an interpretation of the inspected source, not executed CMake output.",
                     "Current tag metadata may differ from historical build metadata; refs and object identities are retained.",
                     "No runtime or scientific setting was changed; this does not explain the separate SIGSEGV."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = inspect(args.repository.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
