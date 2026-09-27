"""Build or restore a bounded scoring archive; raw upstream inputs are excluded."""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tarfile


def identity(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def record(path):
    return dict(path=str(path.absolute()), **identity(path.read_bytes()))


def safe_name(name):
    path = PurePosixPath(name)
    if not path.parts or path.is_absolute() or ".." in path.parts or str(path) != name or "\\" in name:
        raise ValueError("Unsafe archive path")
    return name


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def write_archive(path, payloads):
    with path.open("xb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode="w") as archive:
            for name, data in sorted(payloads.items()):
                member = tarfile.TarInfo(safe_name(name))
                member.size, member.mode, member.mtime = len(data), 0o644, 0
                archive.addfile(member, io.BytesIO(data))


def read_archive(path, sha256):
    data = path.read_bytes()
    if identity(data)["sha256"] != sha256:
        raise ValueError("Archive checksum mismatch")
    payloads, total = {}, 0
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        for member in archive:
            name = safe_name(member.name)
            total += member.size
            if not member.isfile() or name in payloads or member.size < 0 or total > 256 * 1024 * 1024:
                raise ValueError("Unexpected archive member or size")
            payloads[name] = archive.extractfile(member).read()
    manifest = json.loads(payloads["archive_manifest.json"])
    if (manifest["schema_version"] != 1 or manifest["scope"] != "orthobench_scoring_only"
            or manifest["publication_ready"] is not False or manifest["raw_upstream_inputs_included"] is not False):
        raise ValueError("Unexpected archive scope")
    expected = manifest["files"]
    if set(payloads) != set(expected) | {"archive_manifest.json"}:
        raise ValueError("Archive inventory mismatch")
    for name, item in expected.items():
        if identity(payloads[name]) != item:
            raise ValueError("Archive member checksum mismatch")
    return payloads, manifest


def build(repo, score_root, receipt_path, output):
    from benchmark_tools.export_publication_readers import verify
    from benchmark_tools.reproduce_portable_ob_score import verify_files
    if output.exists():
        raise FileExistsError(output)
    receipt = json.loads(receipt_path.read_text())
    if receipt["status"] != "portable_full_orthobench_score_exactly_reproduced" or not receipt["full_score_equal"]:
        raise ValueError("Require the validated full-data scoring receipt")
    for item in [receipt["source"], receipt["reader_manifest"], *receipt["files"]]:
        if record(Path(item["path"])) != item:
            raise ValueError("Changed scoring evidence")
    if {Path(r["path"]).parent for r in receipt["files"]} != {score_root}:
        raise ValueError("Wrong scoring directory")
    reader_manifest = verify(score_root / "readers")
    inputs = json.loads((score_root / "inputs.json").read_text())
    verify_files(score_root, inputs)
    names = ["worker.py", "inputs.json", "data/root_hogs.tsv", "readers/manifest.json",
             *["readers/" + n for n in reader_manifest["files"]]]
    payloads = {name: (score_root / name).read_bytes() for name in names}
    payloads["expected_score.json"] = (score_root / "score.json").read_bytes()
    payloads["requirements.txt"] = (repo / "benchmark_tools/results/publication_reader_requirements_20260927_v2.txt").read_bytes()
    payloads["LICENSE.md"] = (repo / "LICENSE.md").read_bytes()
    payloads["restore.py"] = Path(__file__).read_bytes()
    manifest = dict(schema_version=1, scope="orthobench_scoring_only", publication_ready=False,
        raw_upstream_inputs_included=False, redistribution_cleared=False,
        files={n: identity(p) for n, p in sorted(payloads.items())},
        upstream_revision=inputs["upstream_revision"],
        limitations=["No inference, raw upstream data, wheels, interpreter or OS libraries included.",
                     "Local archive only; no public deposition or redistribution clearance."])
    payloads["archive_manifest.json"] = encoded(manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_archive(output, payloads)
    restored, _ = read_archive(output, record(output)["sha256"])
    if restored != payloads:
        raise ValueError("Archive roundtrip differs")
    return dict(status="portable_ob_scoring_archive_built", archive=record(output),
        receipt=record(receipt_path), files=len(payloads), source=record(Path(__file__)),
        raw_upstream_inputs_included=False, publication_ready=False, redistribution_cleared=False)


def restore(archive, sha256, acquisition, python, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    payloads, manifest = read_archive(archive, sha256)
    inputs = json.loads(payloads["inputs.json"])
    copies = []
    for row in inputs["files"]:
        name = safe_name(row["relative_path"])
        expected = {k: row[k] for k in ("bytes", "sha256")}
        if name in payloads:
            if identity(payloads[name]) != expected:
                raise ValueError("Bundled prediction differs from input manifest")
        else:
            if not name.startswith("data/BENCHMARKS/") or row["role"] not in {"fasta", "reference", "uncertain"}:
                raise ValueError("Unknown external input")
            path = acquisition / name.removeprefix("data/")
            if path.is_symlink() or not path.resolve().is_relative_to(acquisition.resolve()) or identity(path.read_bytes()) != expected:
                raise ValueError("Acquired input differs from archived pin")
            copies.append((path, name, expected))
    output.mkdir(parents=True)
    for name, data in payloads.items():
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    for source, name, expected in copies:
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        if identity(target.read_bytes()) != expected:
            raise ValueError("Acquired input changed during restoration")
    env = dict(HOME="/tmp", PATH="/usr/bin:/bin", LANG="C.UTF-8", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    command = [str(python), "-I", "-B", str(output / "worker.py"), "--worker", str(output),
               "--output", str(output / "restored_score.json")]
    with (output / "restore.log").open("x") as log:
        result = subprocess.run(command, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=300)
    if result.returncode:
        (output / "failure.json").write_bytes(encoded(dict(returncode=result.returncode, command=command, retry=False)))
        raise RuntimeError("Restored scoring failed; no automatic retry")
    observed = json.loads((output / "restored_score.json").read_text())
    expected = json.loads(payloads["expected_score.json"])
    if observed != expected:
        raise ValueError("Restored full score object differs")
    for name, data in payloads.items():
        if (output / name).read_bytes() != data:
            raise ValueError("Archived file changed during restored scoring")
    report = dict(status="archived_full_ob_score_restored_exactly", archive=record(archive),
        upstream_revision=manifest["upstream_revision"], external_files=len(copies),
        command=command, environment=env, score=record(output / "restored_score.json"),
        full_score_equal=True, genes=observed["genes"], groups=observed["groups"],
        refogs=observed["score"]["refogs"], publication_ready=False, inference_rerun=False,
        limitations=["Caller supplies an independently installed reader interpreter and acquired reference data.",
            "Same-host scoring restoration, not fresh inference or cross-platform verification.",
            "No public deposition or redistribution rights clearance."])
    (output / "restoration.json").write_bytes(encoded(report))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="mode", required=True)
    builder = commands.add_parser("build")
    for name in ("repo", "score-root", "receipt", "output"):
        builder.add_argument("--" + name, type=Path, required=True)
    restorer = commands.add_parser("restore")
    for name in ("archive", "acquisition", "python", "output"):
        restorer.add_argument("--" + name, type=Path, required=True)
    restorer.add_argument("--sha256", required=True)
    args = parser.parse_args()
    if args.mode == "build":
        result = build(args.repo.resolve(), args.score_root.resolve(), args.receipt.resolve(), args.output.absolute())
    else:
        result = restore(args.archive.resolve(), args.sha256, args.acquisition.resolve(), args.python.absolute(), args.output.absolute())
    print(json.dumps(result, indent=2, sort_keys=True))
