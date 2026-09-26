"""Bind a setup-overlay installation to frozen scientific Git blobs."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.parse import unquote, urlsplit
import zipfile

from benchmark_tools.prepare_frozen_build_overlay import BUILD_REVISION
from benchmark_tools.verify_frozen_source_archive import REVISION, git_files
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def scientific_members(wheel, files):
    expected = {name: item for name, item in files.items() if name.startswith("orthohmm/")
                and not name.startswith("orthohmm/search/_experimental/")}
    native = {"orthohmm/search/csrc/" + name for name in ("hmm_viterbi.so", "kmer_prefilter.so", "pair_align.so")}
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate wheel member")
        package = {name for name in names if name.startswith("orthohmm/") and not name.endswith("/")}
        if package != set(expected) | native:
            raise ValueError("Unexpected frozen wheel package inventory")
        rows = []
        for name, item in sorted(expected.items()):
            content = archive.read(name)
            if content != item["content"]:
                raise ValueError("Wheel scientific source differs: " + name)
            rows.append(dict(path=name, git_blob=item["git_blob"], bytes=len(content),
                             sha256=hashlib.sha256(content).hexdigest()))
    return rows


def local_install_wheels(report, directory):
    rows = []
    for item in report["install"]:
        url = urlsplit(item["download_info"]["url"])
        if url.scheme != "file" or url.netloc:
            raise ValueError("Require local wheel installation")
        path = Path(unquote(url.path)).resolve()
        if path.parent != directory.resolve() or path.suffix != ".whl":
            raise ValueError("Wheel outside dedicated wheelhouse")
        identity = record(path)
        if identity["sha256"] != item["download_info"]["archive_info"]["hashes"]["sha256"]:
            raise ValueError("Installed wheel hash differs")
        rows.append(dict(name=item["metadata"]["name"], version=item["metadata"]["version"], wheel=identity))
    if (len(rows) != len({r["wheel"]["path"] for r in rows})
            or {Path(r["wheel"]["path"]) for r in rows} != {p.resolve() for p in directory.iterdir()}):
        raise ValueError("Wheelhouse inventory differs from installation")
    return rows


def audit(repo, directory):
    paths = [directory / name for name in ("staging.json", "install_report_v2.json", "smoke/verification.json")]
    inputs = [record(path) for path in paths]
    stage, install, smoke = [json.loads(p.read_text()) for p in paths]
    if stage["scientific_revision"] != REVISION or stage["build_revision"] != BUILD_REVISION:
        raise ValueError("Wrong scientific/build revisions")
    files = git_files(repo)
    overlay = subprocess.check_output(["git", "-C", str(repo), "show", BUILD_REVISION + ":setup.py"])
    if {r["relative_path"] for r in stage["files"]} != set(files) or len(stage["files"]) != len(files):
        raise ValueError("Staged inventory differs")
    for row in stage["files"]:
        name = row["relative_path"]
        expected = overlay if name == "setup.py" else files[name]["content"]
        if (row["setup_overlay"] != (name == "setup.py")
                or row["original_git_blob"] != files[name]["git_blob"]
                or row["original_sha256"] != hashlib.sha256(files[name]["content"]).hexdigest()
                or (directory / "source" / name).read_bytes() != expected):
            raise ValueError("Staged sources no longer match Git")
        inputs.append(row["staged"])
    wheel = Path(smoke["wheel"]["path"])
    check(smoke["wheel"])
    members = scientific_members(wheel, files)
    wheels = local_install_wheels(install, directory / "wheels")
    if smoke["status"] != "same_host_cpu_wheel_smoke_verified" or len(wheels) != 11:
        raise ValueError("Incomplete installed smoke or wheelhouse")
    for row in members:
        installed = Path(smoke["runtime"]["module"]).parent.parent / row["path"]
        if record(installed)["sha256"] != row["sha256"]:
            raise ValueError("Installed scientific bytes differ")
    inputs.extend([smoke["wheel"], *smoke["installed"], *smoke["inputs"], *smoke["original_inputs"],
        *[r["partition"] for r in smoke["runs"]], *[r["log"] for r in smoke["runs"]],
        *[r["wheel"] for r in wheels], stage["historical_setup"], stage["source"], smoke["source"],
        record(directory / "build.log"), record(directory / "install.log"), record(directory / "install_v2.log"),
        record(repo / "benchmark_tools/results/publication_frozen_overlay_requirements_20260926.txt"), record(__file__)])
    check(stage["historical_setup"])
    if Path(stage["historical_setup"]["path"]).read_bytes() != files["setup.py"]["content"]:
        raise ValueError("Historical setup changed")
    for item in inputs:
        check(item)
    return dict(status="frozen_scientific_sources_installed_with_setup_overlay", scientific_revision=REVISION,
        build_revision=BUILD_REVISION, scientific_members=members, wheels=wheels, runs=smoke["runs"],
        checked_records=inputs, original_runtime_reproduced=False, publication_ready=False,
        limitations=["Setup-only overlay and new compilation/dependencies, not the original scientific executable.",
            "Same-host standard/high-sensitivity fixture, not full phylogeny or benchmark equivalence.",
            "Experimental reference sources are excluded from the wheel, but retained in the source archive.",
            "No public release, complete redistribution review or cross-host portability claim."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve(), args.directory.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
