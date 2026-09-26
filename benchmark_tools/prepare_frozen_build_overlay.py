"""Stage frozen scientific sources with an explicitly versioned setup-only overlay."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess

from benchmark_tools.verify_frozen_source_archive import REVISION, git_files
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

BUILD_REVISION = "6fd6df19daba83ec6467b917988f99e27a95be14"


def prepare(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    files = git_files(repo)
    if "setup.py" not in files or not any(name.startswith("orthohmm/") for name in files):
        raise ValueError("Missing package or original setup")
    if any(PurePosixPath(name).is_absolute() or ".." in PurePosixPath(name).parts for name in files):
        raise ValueError("Unsafe source path")
    setup = subprocess.check_output(["git", "-C", str(repo), "show", BUILD_REVISION + ":setup.py"])
    setup_blob = subprocess.check_output(["git", "-C", str(repo), "rev-parse", BUILD_REVISION + ":setup.py"], text=True).strip()
    compile(setup, "setup.py", "exec")
    output.mkdir(parents=True)
    source = output / "source"
    rows = []
    for name, item in sorted(files.items()):
        path = source / name
        path.parent.mkdir(parents=True, exist_ok=True)
        content = setup if name == "setup.py" else item["content"]
        path.write_bytes(content)
        path.chmod(item["mode"])
        rows.append(dict(relative_path=name, original_git_blob=item["git_blob"],
            original_sha256=hashlib.sha256(item["content"]).hexdigest(),
            setup_overlay=name == "setup.py", staged=record(path)))
    original = output / "historical_setup.py"
    original.write_bytes(files["setup.py"]["content"])
    for row in rows:
        check(row["staged"])
        if not row["setup_overlay"] and row["staged"]["sha256"] != row["original_sha256"]:
            raise ValueError("Scientific source changed during staging")
    report = dict(status="frozen_scientific_source_with_setup_overlay_staged",
        scientific_revision=REVISION, build_revision=BUILD_REVISION, setup_blob=setup_blob,
        files=rows, historical_setup=record(original), source=record(__file__),
        publication_ready=False, build_executed=False,
        limitations=["Only setup.py is replaced; scientific sources and version.py remain frozen.",
            "A future wheel is a new compiled artifact, not the historical scientific executable.",
            "No installation, inference, numerical-equivalence or portability claim follows from staging."])
    with (output / "staging.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.repo.resolve(), args.output.absolute())
