"""Recompute the full retained OrthoBench score from independently copied inputs."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

READBACK_SHA = "339dc8e3f72bbfc657a86856f7c5aaa9c29fad1a49537889f1a09e3c3c5acc35"
SCORER_SHA = "0e24c302a9d12ce8c82a4111ff0321c87303469b3c866c317dc2629ad66f9746"


def identity(path):
    data = path.read_bytes()
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def validate_reference_ids(references, uncertain, universe):
    # Frozen low-certainty sets retain blank lines; these are not gene IDs.
    ids = set().union(*references.values(), *uncertain.values()) - {""}
    if not ids <= universe:
        raise ValueError("Reference contains unknown nonempty gene IDs")


def verify_files(root, manifest):
    names = set()
    for row in manifest["files"]:
        name = Path(row["relative_path"])
        if name.is_absolute() or ".." in name.parts or str(name) in names:
            raise ValueError("Unsafe or duplicate staged input")
        names.add(str(name))
        path = root / name
        if (path.is_symlink() or not path.resolve().is_relative_to(root.resolve())
                or identity(path) != {k: row[k] for k in ("bytes", "sha256")}):
            raise ValueError("Changed staged input")
    actual = {str(p.relative_to(root)) for p in (root / "data").rglob("*") if p.is_file()}
    if actual != names:
        raise ValueError("Staged data inventory differs")


def worker(root, output):
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads((root / "inputs.json").read_text())
    verify_files(root, manifest)
    sys.path.insert(0, str(root / "readers"))
    from benchmark_tools.audit_installed_orthobench import read_root_hogs
    from benchmark_tools.run_installed_orthobench import fasta_ids
    from benchmark_tools.score_orthobench_partition import score_partition
    source = root / "readers/benchmark_tools/score_orthobench_partition.py"
    if identity(source)["sha256"] != SCORER_SHA:
        raise ValueError("Changed frozen scorer")
    refs, uncertain, fasta, predictions = {}, {}, [], []
    for row in manifest["files"]:
        path = root / row["relative_path"]
        role = row["role"]
        if role in ("reference", "uncertain"):
            target = refs if role == "reference" else uncertain
            if path.name in target:
                raise ValueError("Duplicate reference name")
            target[path.name] = set(path.read_text().splitlines())
        elif role == "fasta":
            fasta.append(path)
        elif role == "predictions":
            predictions.append(path)
        else:
            raise ValueError("Unknown staged input role")
    if (len(refs), len(uncertain), len(fasta), len(predictions)) != (70, 11, 12, 1):
        raise ValueError("Not the complete OrthoBench input inventory")
    universe = fasta_ids(sorted(fasta))
    if len(universe) != 251378:
        raise ValueError("Wrong input gene universe")
    validate_reference_ids(refs, uncertain, universe)
    groups = read_root_hogs(predictions[0], universe)
    score = score_partition(groups, refs, uncertain)
    verify_files(root, manifest)
    save(output, dict(status="full_orthobench_score_recomputed", genes=len(universe), groups=len(groups),
                      score=score, scorer=identity(source), inputs=identity(root / "inputs.json")))


def reproduce(repo, acquisition, readers, python, output):
    from benchmark_tools.export_publication_readers import verify
    from benchmark_tools.verify_orthobench_acquisition import verify as verify_acquisition
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    upstream = verify_acquisition(acquisition, repo / "benchmark_tools/results")
    reader_manifest = verify(readers)
    readback = repo / "benchmarks/work/publication_full_recovery_orthobench_20260926/readback/result.json"
    if identity(readback)["sha256"] != READBACK_SHA:
        raise ValueError("Changed admitted full-data score")
    retained = json.loads(readback.read_text())
    prediction = retained["partitions"]["current"]
    check(prediction)
    output.mkdir(parents=True)
    shutil.copytree(readers, output / "readers")
    shutil.copyfile(__file__, output / "worker.py")
    rows = []
    for row in upstream["files"]:
        name = row["relative_path"]
        if not name.startswith(("BENCHMARKS/Input/", "BENCHMARKS/RefOGs/")):
            continue
        role = "fasta" if name.startswith("BENCHMARKS/Input/") else (
            "uncertain" if "/low_certainty_assignments/" in name else "reference")
        target = output / "data" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(acquisition / name, target)
        if identity(target) != {k: row[k] for k in ("bytes", "sha256")}:
            raise ValueError("Copy differs from acquired upstream input")
        rows.append(dict(relative_path=str(target.relative_to(output)), role=role, **identity(target)))
    target = output / "data/root_hogs.tsv"
    shutil.copyfile(prediction["path"], target)
    if identity(target) != {k: prediction[k] for k in ("bytes", "sha256")}:
        raise ValueError("Prediction copy differs")
    rows.append(dict(relative_path="data/root_hogs.tsv", role="predictions", **identity(target)))
    save(output / "inputs.json", dict(files=rows, upstream_revision=upstream["upstream_revision"]))
    trace = output / "worker.strace"
    command = ["/usr/bin/strace", "-f", "-s", "4096", "-e", "trace=%file", "-o", str(trace),
               str(python), "-I", "-B", str(output / "worker.py"), "--worker", str(output),
               "--output", str(output / "score.json")]
    environment = dict(HOME="/tmp", PATH="/usr/bin:/bin", LANG="C.UTF-8",
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with (output / "worker.log").open("x") as log:
        result = subprocess.run(command, cwd=output, env=environment, stdout=log,
                                stderr=subprocess.STDOUT, timeout=300)
    if result.returncode:
        save(output / "failure.json", dict(returncode=result.returncode, command=command, retry=False))
        raise RuntimeError("Portable scoring failed; preserve outputs, no automatic retry")
    observed = json.loads((output / "score.json").read_text())
    if observed["score"] != retained["scores"]["current"]:
        raise ValueError("Full score object differs from retained admission")
    if verify(output / "readers") != reader_manifest:
        raise ValueError("Reader export changed")
    verify_files(output, json.loads((output / "inputs.json").read_text()))
    forbidden = [str(repo), "/home/bizon/anaconda3/lib/python3.10/site-packages"]
    found = [p for p in forbidden if p in trace.read_text()]
    if found:
        raise ValueError("Unexpected original path in worker file trace: " + repr(found))
    return dict(status="portable_full_orthobench_score_exactly_reproduced", genes=observed["genes"],
        groups=observed["groups"], refogs=observed["score"]["refogs"], full_score_equal=True,
        scores={k: v for k, v in observed["score"].items() if k != "refog_records"},
        retained=record(readback), prediction=prediction, source=record(__file__),
        reader_manifest=record(output / "readers/manifest.json"), command=command,
        environment=environment, files=[record(output / n) for n in ("worker.py", "inputs.json", "score.json", "worker.log", "worker.strace")],
        forbidden_literal_prefixes=forbidden, forbidden_prefixes_found=found,
        inference_rerun=False, redistribution_cleared=False, publication_ready=False,
        limitations=["Same-host full-data scoring only; no new inference or independent validation.",
            "Original scheduler admission is retained evidence, not rerun by the portable worker.",
            "File trace prefix checks are not a sandbox or OS isolation proof.",
            "Raw benchmark inputs are private acquisition copies, not a redistributable archive."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    for name in ("repo", "acquisition", "readers", "python", "report"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker.resolve(), args.output.absolute())
    else:
        if any(getattr(args, n) is None for n in ("repo", "acquisition", "readers", "python", "report")):
            parser.error("Require --repo, --acquisition, --readers, --python and --report")
        if args.report.exists():
            raise FileExistsError(args.report)
        result = reproduce(*(getattr(args, n).absolute() for n in ("repo", "acquisition", "readers", "python", "output")))
        save(args.report, result)
        print(result["status"])
